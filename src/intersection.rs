use crate::ray::Ray;
use crate::sphere::Sphere;
use crate::tuple::Tuple;


#[macro_export]
macro_rules! agg_intersection{
    ($($intersection:expr),+ $(,)?)=> {
        { // variadic amount of intersections, as well as comma seperator.
            let agg_int = vec![$($intersection:expr),+$(,)?];
            agg_int.into_iter().
            .min_by(|a,b| a.m_time.partial_cmp(&b,m_time).unwrap_or(std::cmp::Ordering::Equal))
        }
    };
}

#[derive (Clone)]
pub struct Intersection<T>{
    m_time : f32,
    m_object : T,
}

impl<T> Intersection<T>{
    pub fn new(t: f32, o : T) -> Self{
        Intersection{
            m_time : t,
            m_object: o,
        }
    }
}

pub fn construe_agg_structs<T: Clone>(_times : Vec<f32>, object : T) -> Vec<Intersection<T>>{
    let mut intersect_agg : Vec<Intersection<T>> = vec![];
    for time in _times{
        if time != -1.{
            let _a = Intersection::new(time, object.clone());
            intersect_agg.push(_a);
        }
    }
    intersect_agg
}

mod test{
    use crate::ray::intersect_unwrap;
    use super::*;
    #[test]
    pub fn agg_intersection(){
        let _r = Ray::new(Tuple::new_point(0.,0.,-5.), Tuple::new_vector(0.,0.,1.));
        let _s = Sphere::new(Tuple::new_point(0.,0.,0.),1. );
        let _xs = intersect_unwrap(&_s, &_r); // so from here
        let _b = construe_agg_structs(_xs, _s.clone()); // to here there has to be a better way?
        assert_eq!(_b[0].m_object, _s);
        assert_eq!(_b[1].m_object, _s);
    }
}
