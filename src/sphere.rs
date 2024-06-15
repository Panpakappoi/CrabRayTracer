use crate::tuple::Tuple;
pub struct Sphere{
    origin: Tuple,
    radius : f32,
}

impl Sphere{
    pub fn new(o: Tuple, r: f32)-> Self{
        Sphere{
            origin :o,
            radius : r,
        }
    }
    pub fn get_origin(&self) -> Tuple{
        self.origin
    }
    pub fn get_radius(&self) -> f32{
        self.radius
    }
}
