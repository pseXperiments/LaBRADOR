use ring_math::custom_ring;
use ring_math::polynomial_ring;
use ring_math::FieldElement;
use ring_math::Polynomial;
use ring_math::PolynomialRingElement;

// creates a scalar ring struct DilithiumRingElement
custom_ring!(DilithiumRingElement, 8380417, "dilithium_23_bit");

// creates a polynomial ring struct
polynomial_ring!(
    DilithiumPolynomialRingElement,
    DilithiumRingElement,
    {
        // creating the ring modulus polynomial
        // here we use x^64 + 1
        let mut p = Polynomial::identity();
        p.term(&DilithiumRingElement::one(), 64);
        p
    },
    "dilithium_x64+1"
);

#[cfg(test)]
mod test {
    use ring_math::Matrix2D;
    use ring_math::Vector;

    use super::*;

    fn sum<T: FieldElement>(v: Vector<T>) -> T {
        v.iter().fold(T::zero(), |acc, x| acc + x.clone())
    }

    #[test]
    fn generate_constraints() {
        // in the paper n = row count and r = column count ?
        let n = 8;
        let r = 8;
        // TODO: clean up rand/sample function naming
        let phi =
            Matrix2D::<DilithiumPolynomialRingElement>::rand_uniform(n, r, &mut rand::thread_rng());
        // let a = DilithiumPolynomialRingElement::sample_rand(&mut rand::thread_rng());
        let a =
            Matrix2D::<DilithiumPolynomialRingElement>::rand_uniform(n, r, &mut rand::thread_rng());
        // unbounded solution for now
        let s =
            Matrix2D::<DilithiumPolynomialRingElement>::rand_uniform(n, r, &mut rand::thread_rng());
        // inner product of si, sj
        let iv = 2;
        let jv = r;
        let mut term1 = DilithiumPolynomialRingElement::zero();
        for i in 0..iv {
            for j in 0..jv {
                let s_inner_prod = sum(s.row(i) * s.row(j));
                // TODO: clean up this matrix syntax
                term1 += s_inner_prod * a.row(i)[j].clone();
            }
        }
        let mut term2 = DilithiumPolynomialRingElement::zero();
        for i in 0..r {
            term2 += sum(phi.row(i).clone() * s.row(i).clone());
        }
        // solve for a b, then check equation to 0
        let b = -(term1.clone() + term2.clone());
        println!("term1: {term1}");
        println!("term2: {term2}");
        println!("b: {b}");
        assert_eq!(term1 + term2 + b, DilithiumPolynomialRingElement::zero());
    }
}
