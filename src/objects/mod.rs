use crate::{intersect::Intersection, Ray};

pub mod sphere;

pub trait Object {
    fn intersect(&self, r: &Ray, tmin: f32, tmax: f32) -> Option<Intersection>;
}
