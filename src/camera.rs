use std::{f32::consts::PI, fs::File, io::Write};

use fastrand::Rng;
use glam::Vec3A;

use crate::{intersect::Intersection, material::ReflT, scene::Scene};

pub struct Ray {
    pub o: Vec3A,
    pub d: Vec3A,
}

impl Ray {
    pub fn at(&self, t: f32) -> Vec3A {
        return self.o + self.d * t;
    }
}

pub struct Camera {
    pub position: Vec3A,
    pub point_to: Vec3A,
}

impl Camera {
    pub fn capture(&self, width: usize, height: usize, samps: usize, scene: &Scene) {
        let cx = Vec3A::new(width as f32 * 0.5135 / height as f32, 0., 0.);
        let cy = cx.cross(self.position).normalize() * 0.5135;
        let mut c = vec![Vec3A::ZERO; width * height];
        let mut rng = Rng::new();
        for y in 0..height {
            print!(
                "\rRendering ({} spp) {}",
                samps * 4,
                100. * y as f32 / (height as f32 - 1.)
            );
            for x in 0..width {
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
                            let d = cx
                                * (((sx as f32 + 0.5 + dx) / 2. + x as f32) / width as f32 - 0.5)
                                + cy * (((sy as f32 + 0.5 + dy) / 2. + y as f32) / height as f32
                                    - 0.5)
                                + self.point_to;
                            let ray = Ray {
                                o: self.position + d * 140.,
                                d: d.normalize(),
                            };
                            r += radiance(&ray, &scene, 0, &mut rng) / samps as f32;
                        }
                        let i = y * width + x;
                        c[i] = c[i] + r.clamp(Vec3A::ZERO, Vec3A::ONE) / 4.;
                    }
                }
            }
        }

        let mut file = File::create("result/image.ppm").unwrap();
        file.write(format!("P3\n{} {}\n{}\n", width, height, 255).as_bytes())
            .unwrap();
        for color in c.into_iter() {
            let (x, y, z) = (color.x, color.y, color.z);
            file.write(format!("{} {} {} ", to_int(x), to_int(y), to_int(z)).as_bytes())
                .unwrap();
        }
    }
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

fn to_int(x: f32) -> i32 {
    return (x.clamp(0., 1.).powf(1. / 2.2) * 255. + 0.5) as i32;
}
