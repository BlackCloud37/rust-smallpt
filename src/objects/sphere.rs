use glam::Vec3A;

use crate::{camera::Ray, ReflT};

use super::{Intersection, Object};
#[derive(Debug)]
pub struct Sphere {
    rad: f32,
    p: Vec3A,
    e: Vec3A,
    c: Vec3A,
    refl: ReflT,
}

impl Sphere {
    pub fn new(rad: f32, p: Vec3A, e: Vec3A, c: Vec3A, refl: ReflT) -> Self {
        Self { rad, p, e, c, refl }
    }
}

impl Object for Sphere {
    fn intersect(&self, r: &Ray, tmin: f32, tmax: f32) -> Option<Intersection> {
        let op = self.p - r.o;
        let b = op.dot(r.d);
        let mut det = b * b - op.dot(op) + self.rad * self.rad;
        if det < 0. {
            return None;
        } else {
            det = det.sqrt();
        }
        let t = b - det;
        if t >= tmin && t <= tmax {
            let point = r.at(t);
            return Some(Intersection {
                point,
                norm: (point - self.p).normalize(),
                e: self.e,
                reflt: self.refl,
                c: self.c,
                t,
            });
        } else {
            let t = b + det;
            if t >= tmin && t <= tmax {
                let point = r.at(t);
                return Some(Intersection {
                    point,
                    norm: (point - self.p).normalize(),
                    e: self.e,
                    reflt: self.refl,
                    c: self.c,
                    t,
                });
            } else {
                return None;
            }
        }
    }
}
