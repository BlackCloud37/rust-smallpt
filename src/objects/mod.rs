use crate::{camera::Ray, intersect::Intersection};

pub mod sphere;

pub trait Object: Sync + Send {
    fn intersect(&self, r: &Ray, tmin: f32, tmax: f32) -> Option<Intersection>;
}
