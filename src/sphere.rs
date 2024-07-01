use crate::matrix::{identity, Mat};
use crate::tuple::Tuple;
#[derive (Clone, PartialEq, Debug)]
pub struct Sphere{
    origin: Tuple,
    radius : f32,
    transform: Mat,
}


impl Sphere{
    pub fn new(o: Tuple, r: f32, m:Option<Mat>)-> Self{
        Sphere{
            origin :o,
            radius : r,
            transform: m.unwrap_or(identity()).into()
        }
    }
    pub fn get_origin(&self) -> Tuple{
        self.origin
    }
    pub fn get_radius(&self) -> f32{
        self.radius
    }

    pub fn set_transform(&mut self, _m: Mat){
        self.transform = _m;
    }
    pub fn get_transform(&self)->Mat{self.transform.clone()}
}
#[cfg(test)]
mod tests{
    use crate::intersection::construe_agg_structs;
    use crate::matrix::scaling;
    use crate::ray::intersect_unwrap;
    use crate::matrix::translation;
    use crate::ray::Ray;
    use super::*;
    #[test]
    pub fn set_transform_test(){
        let mut _s = Sphere::new(Tuple::new_point(0.,0.,0.),1.,None);
        let _t = translation(2.,3.,4.);
        _s.set_transform(_t);
        let _b = _s.get_transform();
        assert_eq!(_s.transform , translation(2.,3.,4.));
    }
    #[test]
    pub fn intersect_scaled_sphere_with_ray(){
        let mut _r = Ray::new(Tuple::new_point(0., 0., -5.), Tuple::new_vector(0.,0.,1.));
        let mut _s = Sphere::new(Tuple::new_point(1.,1.,1.), 1., None);
        _s.set_transform(scaling(2.,2.,2.));

        let _xs = intersect_unwrap(&_s, &_r);
        let _xs_agg = construe_agg_structs(_xs, _s.clone());
        assert_eq!(_xs_agg.len(), 2);
        assert_eq!(_xs_agg[0].get_time(), 3.);
        assert_eq!(_xs_agg[1].get_time(), 7.);
    }
}