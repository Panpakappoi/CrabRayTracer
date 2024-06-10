use crate::tuple::Tuple;
#[derive(Clone, Copy)]
pub struct Environment{
    pub wind : Tuple,
    pub gravity : Tuple
}
#[derive(Clone, Copy)]
pub struct Projectile{
    pub position : Tuple,
    pub velocity : Tuple
}

impl Projectile{
    pub fn new( x_1:f32, y_1:f32, z_1:f32, x_2:f32, y_2:f32, z_2:f32) -> Self{
        let position = Tuple::new_point(x_1,y_1,z_1);
        let velocity = Tuple::new_vector(x_2,y_2,z_2);
        Projectile{position , velocity}
    }
}

impl Environment{
    pub fn new(x_1 :f32, y_1:f32, z_1:f32, x_2 :f32, y_2:f32, z_2: f32) -> Self{
        let wind = Tuple::new_vector(x_1,y_1,z_1);
        let gravity = Tuple::new_vector(x_2, y_2, z_2);
        Environment{wind, gravity}
    }
}

pub fn tick(env : &Environment, proj : &Projectile) -> Projectile{
    let mut position = proj.position + proj.velocity;
    let mut velocity = proj.velocity + env.wind + env.gravity;
    
    let new_pos = Projectile::new( position.x, position.y, position.z, velocity.x, velocity.y, velocity.z);
    println!("{} {} {} ", new_pos.position.x, new_pos.position.y, new_pos.position.z);
    new_pos
}

