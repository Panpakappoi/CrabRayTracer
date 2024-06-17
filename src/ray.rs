use crate::tuple::Tuple;
use crate::sphere::Sphere;
pub fn intersect_unwrap(sphere : &Sphere, ray:&Ray) -> Vec<f32>{
    intersect(sphere, ray).unwrap_or(vec![-1.])
}

pub fn intersect(sphere : &Sphere, ray:&Ray)-> Option<Vec<f32>>{
    let mut intersections = vec![] as Vec<f32>;
    let s_origin = sphere.get_origin();
    let origin_to_center = Tuple::new_vector(
        ray.m_origin.get_X() - s_origin.get_X(),
        ray.m_origin.get_Y() - s_origin.get_Y(),
        ray.m_origin.get_Z() - s_origin.get_Z(),
    );
    let _a = ray.m_direction.dot_product(ray.m_direction);
    println!("{}", _a);
    let _b = 2.0 * origin_to_center.dot_product(ray.m_direction);
    println!("{}", _b);
    let _c = origin_to_center.dot_product(origin_to_center) - (sphere.get_radius() * sphere.get_radius());
    println!("{}", _c);
    let _disc  = (_b * _b) - (4.0*_a*_c);
    println!("{}", _disc);
    if _disc < 0. {
        None
    } else {
        let _t1 = (-_b - _disc.sqrt()) / (2.0 * _a);
        let _t2 = (-_b + _disc.sqrt()) / (2.0 * _a);
        if _t1 < 0.{
            intersections.push(_t2);
            return Some(intersections);
        }
        if _t2 < 0.{
            intersections.push(_t1);
            return Some(intersections);
        }
        else {
            intersections.push(_t1);
            intersections.push(_t2);
            return Some(intersections);
        }
    }

 }

pub fn aggregate_intersection(){

}
pub struct Ray{
    m_origin: Tuple,
    m_direction: Tuple,
}

impl Ray{
    pub fn new(origin: Tuple, direction: Tuple)-> Self{
        Ray{
            m_origin: origin,
            m_direction : direction
        }
    }
    pub fn position(&self, time : f32) -> Tuple {
        self.m_origin + self.m_direction * time
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn create_ray_test(){
        let _origin = Tuple::new_point(1.,2.,3.);
        let _direction = Tuple::new_vector(4.,5.,6.);
        let _a = Ray::new(
            Tuple::new_point(1.,2.,3.),
            Tuple::new_vector(4.,5.,6.));
        assert_eq!(_a.m_origin, _origin);
        assert_eq!(_a.m_direction, _direction);
    }
    #[test]
    fn compute_point_dist(){
        let _r = Ray::new(
            Tuple::new_point(2.,3.,4.),
            Tuple::new_vector(1.,0.,0.));
        assert_eq!(_r.position(0.), Tuple::new_point(2.,3.,4.));
        assert_eq!(_r.position(1.), Tuple::new_point(3.,3.,4.));
        assert_eq!(_r.position(-1.), Tuple::new_point(1.,3.,4.));
        assert_eq!(_r.position(2.5), Tuple::new_point(4.5,3.,4.));
    }
    #[test]
    fn intersection(){
        let _r = Ray::new(Tuple::new_point(0.,0.,-5.), Tuple::new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.), 1.);
        let _xs = intersect(&_s, &_r);
        let _t = _xs.unwrap();
        assert_eq!(_t[0] , 4.);
        assert_eq!(_t[1], 6.);

    }


}