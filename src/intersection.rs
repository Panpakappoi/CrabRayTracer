use crate::ray::Ray;
use crate::sphere::Sphere;
use crate::tuple::Tuple;

/*
Problem with scalability. Think about many different intersections. On multiple objects. Kind of Just set up as simple test on objects.
Would it make sense to store. I mean... Have to check in 3 dimensions for an intersection?

So it would be ideal if we were able to just hashmap to index to the coordinates of the object. Key are the coordinates.

 */
// This is so jank...
fn calc_hit<T:Clone>(intersections: Vec<Intersection<T>>) -> Option<Intersection<T>>{
    let mut e_hit = <f32>::INFINITY;
    let mut first_intersection :Option<Intersection<T>> = None;
    let mut i = intersections.len();
    for intersect in intersections{
        if intersect.m_time > 0. && intersect.m_time < e_hit{
            e_hit = intersect.m_time;
            first_intersection = Some(intersect);
        }
    }
    if e_hit != <f32>::INFINITY{
        return first_intersection;
    }
        None

}


// #[macro_export]
// macro_rules! agg_intersection{
//     ($($intersection:expr),+ $($intersection:expr)?)=> {
//         { // variadic amount of intersections, as well as comma seperator.
//             let agg_int = vec![$($intersection:expr),+$($intersection:expr)?];
//             agg_int.into_iter()
//             .min_by(|a,b| a.m_time.partial_cmp(&b,m_time).unwrap_or(std::cmp::Ordering::Equal))
//         }
//     };
// }

#[derive (Clone, Debug, PartialEq)]
pub struct Intersection<T>{
    m_time : f32,
    m_object : T,
}


impl<T: Clone> Intersection<T>{
    pub fn new(t: f32, o : T) -> Self{
        Intersection{
            m_time : t,
            m_object: o,
        }
    }
    pub fn get_time(&self) -> f32 {
       self.m_time
    }
    pub fn get_object(&self) -> T{
        self.m_object.clone()
    }
}

pub fn construe_agg_structs<T: Clone>(_times : Vec<f32>, object : T) -> Vec<Intersection<T>>{
    let mut intersect_agg : Vec<Intersection<T>> = vec![];
    for time in _times{
        let _a = Intersection::new(time, object.clone());
        intersect_agg.push(_a);
    }
    intersect_agg
}

mod test{
    use crate::ray::intersect;
    use crate::ray::intersect_unwrap;
    use super::*;
    #[test]
    pub fn agg_intersection_test(){
        let _r = Ray::new(Tuple::new_point(0.,0.,-5.), Tuple::new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1., None);
        let _xs = intersect_unwrap(&_s, &_r); // so from here
        let _b = construe_agg_structs(_xs, _s.clone()); // to here there has to be a better way?
        assert_eq!(_b[0].m_time, 4.);
        assert_eq!(_b[1].m_time, 6.);
    }
    #[test]
    pub fn intersect_sphere_tan_test(){
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1., None);
        let _r = Ray::new(Tuple::new_point(0.,1.,-5.),Tuple::new_vector(0.,0.,1.));
        let _xs = intersect_unwrap(&_s, &_r);
        let _xs_agg = construe_agg_structs(_xs, _s.clone());
        assert_eq!(_xs_agg[0].m_time, 5.);
        assert_eq!(_xs_agg[1].m_time, 5.);
    }
    #[test]
    pub fn intersect_miss_test(){
        let _r = Ray::new(Tuple::new_point(0., 2., -5.), Tuple:: new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1., None);
        let _xs = intersect_unwrap(&_s, &_r);
        let _xs_agg = construe_agg_structs(_xs, _s.clone());
        assert_eq!(_xs_agg.len(), 0);
    }
    #[test]
    pub fn ray_inside_sphere_test(){
        let _r = Ray::new(Tuple::new_point(0.,0.,0.), Tuple::new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1., None);
        let _xs = intersect_unwrap(&_s, &_r);
        let _xs_agg = construe_agg_structs(_xs, _s.clone());
        assert_eq!(_xs_agg[0].m_time, -1.);
        assert_eq!(_xs_agg[1].m_time, 1.);
    }
    #[test]
    pub fn ray_behind_sphere_test(){
        let _r = Ray::new(Tuple::new_point(0., 0., 5.), Tuple::new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1., None);
        let _xs = intersect_unwrap(&_s, &_r);
        let _xs_agg = construe_agg_structs(_xs, _s.clone());
        assert_eq!(_xs_agg[0].m_time, -6.);
        assert_eq!(_xs_agg[1].m_time, -4.);
    }
    #[test]
    pub fn calc_hit_test(){
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.), 1., None);
        let _i1 = Intersection::new(1., &_s);
        let _i2 = Intersection::new(2., &_s);
        let _agg = vec![_i1.clone(), _i2.clone()];
        assert_eq!(calc_hit(_agg).unwrap(), _i1);
    }
    #[test]
    pub fn calc_one_hit_test() {
        let _s = Sphere::new(Tuple::new_point(0., 0., 0.), 1., None);
        let _i1 = Intersection::new(-1., &_s);
        let _i2 = Intersection::new(1., &_s);
        let _agg = vec![_i1.clone(), _i2.clone()];
        assert_eq!(calc_hit(_agg).unwrap(), _i2);
    }
    #[test]
    pub fn calc_hit_no_valid_test(){
        let _s = Sphere::new(Tuple::new_point(0., 0., 0.), 1., None);
        let _i1 = Intersection::new(-1., &_s);
        let _i2 = Intersection::new(-2., &_s);
        let _agg = vec![_i1.clone(), _i2.clone()];
        assert_eq!(calc_hit(_agg), None);
    }
    #[test]
    pub fn calc_hit_multi_test(){
        let _s = Sphere::new(Tuple::new_point(0., 0., 0.), 1., None);
        let _i1 = Intersection::new(5., &_s);
        let _i2 = Intersection::new(7., &_s);
        let _i3 = Intersection::new(-3., &_s);
        let _i4 = Intersection::new(2., &_s);
        let _agg = vec![_i1.clone(), _i2.clone(), _i3.clone(), _i4.clone()];
        assert_eq!(calc_hit(_agg).unwrap(), _i4);
    }

}
