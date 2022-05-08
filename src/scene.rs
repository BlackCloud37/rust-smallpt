use std::f32::INFINITY;

use glam::Vec3A;

use crate::{
    intersect::Intersection,
    objects::{sphere::Sphere, Object},
    Ray,
};

pub struct Scene {
    objs: Vec<Box<dyn Object>>,
}

impl Scene {
    pub fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let mut t = INFINITY;
        let mut h = None;
        for (_, sphere) in self.objs.iter().enumerate() {
            if let Some(hit) = sphere.intersect(r, 0.01, t) {
                t = hit.t;
                h = Some(hit);
            }
        }
        return h;
    }
}

pub fn scene_cornell_box() -> Scene {
    use crate::ReflT::*;
    let spheres: Vec<Sphere> = vec![
        Sphere::new(
            1e5,
            Vec3A::new(1e5 + 1., 40.8, 81.6),
            Vec3A::ZERO,
            Vec3A::new(0.75, 0.25, 0.25),
            DIFF,
        ), //Left
        Sphere::new(
            1e5,
            Vec3A::new(-1e5 + 99., 40.8, 81.6),
            Vec3A::ZERO,
            Vec3A::new(0.25, 0.25, 0.75),
            DIFF,
        ), //Rght
        Sphere::new(
            1e5,
            Vec3A::new(50., 40.8, 1e5),
            Vec3A::ZERO,
            Vec3A::new(0.75, 0.75, 0.75),
            DIFF,
        ), //Back
        Sphere::new(
            1e5,
            Vec3A::new(50., 40.8, -1e5 + 170.),
            Vec3A::ZERO,
            Vec3A::ZERO,
            DIFF,
        ), //Frnt
        Sphere::new(
            1e5,
            Vec3A::new(50., 1e5, 81.6),
            Vec3A::ZERO,
            Vec3A::new(0.75, 0.75, 0.75),
            DIFF,
        ), //Botm
        Sphere::new(
            1e5,
            Vec3A::new(50., -1e5 + 81.6, 81.6),
            Vec3A::ZERO,
            Vec3A::new(0.75, 0.75, 0.75),
            DIFF,
        ), //Top
        Sphere::new(
            16.5,
            Vec3A::new(27., 16.5, 47.),
            Vec3A::ZERO,
            Vec3A::ONE * 0.999,
            SPEC,
        ), //Mirr
        Sphere::new(
            16.5,
            Vec3A::new(73., 16.5, 78.),
            Vec3A::ZERO,
            Vec3A::ONE * 0.999,
            REFR,
        ), //Glas
        Sphere::new(
            600.,
            Vec3A::new(50., 681.6 - 0.27, 81.6),
            Vec3A::new(12., 12., 12.),
            Vec3A::ZERO,
            DIFF,
        ), //Lite
    ];
    let spheres = spheres
        .into_iter()
        .map(|s| Box::new(s) as Box<dyn Object>)
        .collect();
    return Scene { objs: spheres };
}
