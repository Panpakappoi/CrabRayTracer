extern crate RaytracerChallenge;
use RaytracerChallenge::tuple::Tuple;

fn main() {
    let environment = World::new(Tuple::new_vector(0.0_f64, -0.1_f64, 0.0_f64), Tuple::new_vector(-0.01_f64, 0.0_f64,0.0_f64));
    let projectile = Projectile::new(Tuple::new_point(0.0_f64, 20.0_f64, 0.0_f64), Tuple::new_vector(4.0_f64, 5.0_f64, 0.0_f64).normalize());
    let mut current = projectile;
    let mut _tick = 0; 
    while current.position.y > 0.0
    {   
        current = tick(&environment, &current);
        _tick+= 1;
        println!("tick: {}, projectile x_ position {}, projectile y_position {}, projectile z_position{}", _tick, current.position.x, current.position.y,current.position.z);
    }
}


pub struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

impl Projectile {
    pub fn new(position: Tuple, velocity:Tuple) -> Self{
        Projectile {
            position, velocity
        }
    }
}

pub fn tick(env: &World, projectile: &Projectile) -> Projectile {
    Projectile::new(projectile.position + projectile.velocity, projectile.velocity + env.gravity + env.wind)  
}

pub struct World {
    gravity: Tuple,
    wind: Tuple,
}
impl World {
    pub fn new(gravity: Tuple, wind:Tuple) -> Self{
        World{
            gravity, wind
        }
    }
}





