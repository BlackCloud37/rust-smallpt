use fastrand::Rng;
use glam::Vec3A;
use intersect::Intersection;
use material::ReflT;
use scene::{scene_cornell_box, Scene};
use std::{f32::consts::PI, fs::File, io::Write};

mod intersect;
mod material;
mod objects;
mod scene;

#[derive(Debug)]
pub struct Ray {
    pub o: Vec3A,
    pub d: Vec3A,
}

impl Ray {
    pub fn at(&self, t: f32) -> Vec3A {
        return self.o + self.d * t;
    }
}

fn to_int(x: f32) -> i32 {
    return (x.clamp(0., 1.).powf(1. / 2.2) * 255. + 0.5) as i32;
}

fn radiance(r: &Ray, scene: &Scene, depth: usize, rng: &mut Rng) -> Vec3A {
    if let Some(Intersection {
        reflt,
        point: x,
        e,
        c,
        norm: n,
        t: _,
    }) = scene.intersect(r)
    {
        let nl = if n.dot(r.d) < 0. { n } else { -n };
        let mut f = c;
        let p = f.max_element();
        let depth = depth + 1;
        if depth > 10 {
            return e;
        }
        if depth > 5 {
            if rng.f32() < p {
                f /= p;
            } else {
                return e;
            }
        }

        match reflt {
            ReflT::DIFF => {
                let r1 = 2. * PI * rng.f32();
                let r2 = rng.f32();
                let r2s = r2.sqrt();
                let w = nl;
                let u = ((if w.x.abs() > 0.1 {
                    Vec3A::new(0., 1., 0.)
                } else {
                    Vec3A::new(1., 0., 0.)
                })
                .cross(w))
                .normalize();
                let v = w.cross(u);
                let d =
                    (u * r1.cos() * r2s + v * r1.sin() * r2s + w * (1. - r2).sqrt()).normalize();
                let newr = Ray { o: x + 0.1 * d, d };
                return f.mul_add(radiance(&newr, scene, depth, rng), e);
            }
            ReflT::SPEC => {
                let d = r.d - n * 2. * n.dot(r.d);
                let newr = Ray { o: x, d };
                return f.mul_add(radiance(&newr, scene, depth, rng), e);
            }
            ReflT::REFR => {
                let refld = r.d - n * 2. * n.dot(r.d);
                let reflRay = Ray { o: x, d: refld };
                let into = n.dot(nl) > 0.;
                let nc = 1.;
                let nt = 1.5;
                let nnt = if into { nc / nt } else { nt / nc };
                let ddn = r.d.dot(nl);
                let cos2t = 1. - nnt * nnt * (1. - ddn * ddn);
                if cos2t < 0. {
                    return f.mul_add(radiance(&reflRay, scene, depth, rng), e);
                }
                let tdir = (r.d * nnt
                    - n * (if into { 1. } else { -1. }) * (ddn * nnt + cos2t.sqrt()))
                .normalize();
                let a = nt - nc;
                let b = nt + nc;
                let R0 = a * a / b * b;
                let c = 1. - if into { -ddn } else { tdir.dot(n) };
                let Re = R0 + (1. - R0) * c * c * c * c * c;
                let Tr = 1. - Re;
                let P = 0.25 + 0.5 * Re;
                let RP = Re / P;
                let TP = Tr / (1. - P);
                let tray = Ray { o: x, d: tdir };
                return f.mul_add(
                    if depth > 2 {
                        if rng.f32() < P {
                            radiance(&reflRay, scene, depth, rng) * RP
                        } else {
                            radiance(&tray, scene, depth, rng) * TP
                        }
                    } else {
                        radiance(&reflRay, scene, depth, rng) * Re
                            + radiance(&tray, scene, depth, rng) * Tr
                    },
                    e,
                );
            }
        }
    }
    Vec3A::ZERO
}

fn main() {
    let scene = scene_cornell_box();

    let w = 1024 / 2;
    let h = 768 / 2;
    let samps = 4;
    let cam = Ray {
        o: Vec3A::new(50., 52., 295.6),
        d: Vec3A::new(0., -0.042612, -1.).normalize(),
    };
    let cx = Vec3A::new(w as f32 * 0.5135 / h as f32, 0., 0.);
    let cy = cx.cross(cam.d).normalize() * 0.5135;
    let mut c = vec![Vec3A::ZERO; w * h];
    let mut rng = Rng::new();
    for y in 0..h {
        print!(
            "\rRendering ({} spp) {}",
            samps * 4,
            100. * y as f32 / (h as f32 - 1.)
        );
        for x in 0..w {
            for sy in 0..2 {
                for sx in 0..2 {
                    let mut r = Vec3A::ZERO;
                    for _ in 0..samps {
                        let r1 = 2. * rng.f32();
                        let dx = if r1 < 1. {
                            r1.sqrt() - 1.
                        } else {
                            1. - (2. - r1).sqrt()
                        };
                        let r2 = 2. * rng.f32();
                        let dy = if r2 < 1. {
                            r2.sqrt() - 1.
                        } else {
                            1. - (2. - r2).sqrt()
                        };
                        let d = cx * (((sx as f32 + 0.5 + dx) / 2. + x as f32) / w as f32 - 0.5)
                            + cy * (((sy as f32 + 0.5 + dy) / 2. + y as f32) / h as f32 - 0.5)
                            + cam.d;
                        let ray = Ray {
                            o: cam.o + d * 140.,
                            d: d.normalize(),
                        };
                        r += radiance(&ray, &scene, 0, &mut rng) / samps as f32;
                    }
                    let i = (h - y - 1) * w + x;
                    c[i] = c[i] + r.clamp(Vec3A::ZERO, Vec3A::ONE) / 4.;
                }
            }
        }
    }
    let mut file = File::create("result/image.ppm").unwrap();
    file.write(format!("P3\n{} {}\n{}\n", w, h, 255).as_bytes())
        .unwrap();
    for i in 0..w * h {
        file.write(format!("{} {} {} ", to_int(c[i].x), to_int(c[i].y), to_int(c[i].z)).as_bytes())
            .unwrap();
    }
}
