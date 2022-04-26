use glam::Vec3A;

use crate::ReflT;

pub struct Intersection {
    pub t: f32,
    pub point: Vec3A,
    pub norm: Vec3A,

    // Mat
    pub reflt: ReflT,
    pub e: Vec3A,
    pub c: Vec3A,
}
