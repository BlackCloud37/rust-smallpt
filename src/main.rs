use std::sync::Arc;

use camera::Camera;
use glam::Vec3A;
use material::ReflT;
use scene::scene_cornell_box;

mod camera;
mod intersect;
mod material;
mod objects;
mod scene;

fn main() {
    let scene = Arc::new(scene_cornell_box());
    let camera = Camera {
        position: Vec3A::new(50., 52., 295.6),
        point_to: Vec3A::new(0., -0.042612, -1.).normalize(),
    };
    camera.capture(1024 / 2, 768 / 2, 100, scene);
}
