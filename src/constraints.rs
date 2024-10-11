use ring_math::polynomial_ring;
use ring_math::Polynomial;
use ring_math::PolynomialRingElement;
use scalarff::scalar_ring;
use scalarff::FieldElement;

// creates a scalar ring struct DilithiumRingElement
scalar_ring!(DilithiumRingElement, 8380417, "dilithium_23_bit");

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

    use super::*;

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
        let a = DilithiumPolynomialRingElement::sample_uniform(&mut rand::thread_rng());
        let mut a_automorphism = DilithiumPolynomialRingElement::zero();
        for (i, coef) in a.coef().iter().enumerate() {
            if i == 0 {
                // constant term is unchanged
                // x^-1^0 = 1
                a_automorphism.0.term(coef, 0);
                continue;
            }
            let j = degree - i;
            a_automorphism.0.term(&(-(*coef)), j);
        }

        let b = DilithiumPolynomialRingElement::sample_uniform(&mut rand::thread_rng());

        let ct = (a_automorphism * b.clone()).0.constant_term();

        // coefficient dot product of a and b coefficients
        let coef_dot_prod = a.coef().dot_product(b.coef());
        assert_eq!(coef_dot_prod, ct);
    }

    #[test]
    fn generate_constraints() {
        // in the paper n = row count and r = column count ?
        let n = 8;
        let r = 8;
        let phi = Matrix2D::<DilithiumPolynomialRingElement>::sample_uniform(
            n,
            r,
            &mut rand::thread_rng(),
        );
        // let a = DilithiumPolynomialRingElement::sample_uniform(&mut rand::thread_rng());
        let a = Matrix2D::<DilithiumPolynomialRingElement>::sample_uniform(
            n,
            r,
            &mut rand::thread_rng(),
        );
        // unbounded solution for now
        let s = Matrix2D::<DilithiumPolynomialRingElement>::sample_uniform(
            n,
            r,
            &mut rand::thread_rng(),
        );
        // inner product of si, sj
        let iv = 2;
        let jv = r;
        let mut term1 = DilithiumPolynomialRingElement::zero();
        for i in 0..iv {
            for j in 0..jv {
                let s_inner_prod = s.row(i).dot_product(s.row(j));
                // TODO: clean up this matrix syntax
                term1 += s_inner_prod * a.row(i)[j].clone();
            }
        }
        let mut term2 = DilithiumPolynomialRingElement::zero();
        for i in 0..r {
            term2 += phi.row(i).dot_product(s.row(i));
        }
        // solve for a b, then check equation to 0
        let b = -(term1.clone() + term2.clone());
        println!("term1: {term1}");
        println!("term2: {term2}");
        println!("b: {b}");
        assert_eq!(term1 + term2 + b, DilithiumPolynomialRingElement::zero());
    }
}
