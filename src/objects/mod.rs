use glam::Vec3A;

use crate::{Ray, ReflT};

pub mod sphere;

pub struct HitPoint {
    pub t: f32,
    // Material
    pub reflt: ReflT,
    pub p: Vec3A,
    pub e: Vec3A,
    pub c: Vec3A,
}

pub trait Object {
    fn intersect(&self, r: &Ray, tmin: f32, tmax: f32) -> Option<HitPoint>;
}
