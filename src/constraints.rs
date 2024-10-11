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
    fn conj_automorphism() {
        // LaBRADOR page 6, conjugation automorphism
        //
        // polynomial implementation does not allow negative exponents
        // we'll manually find a congruent value for x^-1
        // for now we'll assert the ring modulus is x^n + 1
        // with n >= 2
        let m = DilithiumPolynomialRingElement::modulus();
        let degree = m.degree();
        assert!(degree >= 2);
        for (i, coef) in m.coefficients.iter().enumerate() {
            match i {
                0 => assert_eq!(coef, &DilithiumRingElement::one()),
                _ if i == degree => assert_eq!(coef, &DilithiumRingElement::one()),
                _ => assert_eq!(coef, &DilithiumRingElement::zero()),
            }
        }
        // we'll determine a polynomial with positive exponents
        // that is congruent to x^-1
        // in Z_q[x]/(x^n + 1)
        // -x^(n-1) = x^-1
        // x^-z = -x^(n-z)
        let a = DilithiumPolynomialRingElement::sample_rand(&mut rand::thread_rng());
        let mut out_coefs = vec![DilithiumRingElement::zero(); degree];
        for (i, coef) in a.polynomial().coefficients.iter().enumerate() {
            if i == 0 {
                // constant term is unchanged
                // x^-1^0 = 1
                out_coefs[0] = coef.clone();
                continue;
            }
            let j = degree - i;
            out_coefs[j] = -coef.clone();
        }
        let a_automorphism = DilithiumPolynomialRingElement(Polynomial {
            coefficients: out_coefs,
        });

        let b = DilithiumPolynomialRingElement::sample_rand(&mut rand::thread_rng());

        let ct = (a_automorphism * b.clone()).polynomial().coefficients[0];

        // coefficient dot product of a and b coefficients
        let coef_dot_prod = {
            let mut out = DilithiumRingElement::zero();
            let zero = DilithiumRingElement::zero();
            for i in 0..(usize::max(a.polynomial().degree(), b.polynomial().degree()) + 1) {
                out += a.polynomial().coefficients.get(i).unwrap_or(&zero).clone()
                    * b.polynomial().coefficients.get(i).unwrap_or(&zero).clone();
            }
            out
        };
        assert_eq!(coef_dot_prod, ct);
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
