

/** 
* Evaluates the cumulative distribution function (CDF) for a Student's t distribution with degrees of freedom `v` at a value `t`.
*
* @params {number} t - value for the t test
* @params {PositiveNumber} v - degree of freedom
* @returns {Probability} evaluated CDF
*/
function TdistributionCDF(t, v){
    if (isNaN(t) || isNaN(v) || v <= 0.0) {
        return NaN;
    }
    if (t === 0.0) {
        return 0.5;
    }

    // return 0.5 + 0.5 * Math.sign(t) * (1 - incompbeta(0.5*v, 0.5, v/(v+t*t)));
    let x = v / (v + t * t);
    
    // Handle large x cases to prevent NaN
    if (x < 1e-15) {
        return 1.0; // Approximate result for very large t-values
    }

    let beta_val = incompbeta(0.5 * v, 0.5, x);

    if (isNaN(beta_val)) {
        console.log("Error: incompbeta returned NaN for", x, v);
        return NaN;
    }

    return 0.5 + 0.5 * Math.sign(t) * (1 - beta_val);
}


/** 
* incompbeta(a,b,x) evaluates incomplete beta function, here a, b > 0 and 0 <= x <= 1. This function requires contfractbeta(a,b,x, ITMAX = 200) 
*  code translated from  https://malishoaib.wordpress.com/2014/04/15/the-beautiful-beta-functions-in-raw-python/
*
* @params
* @params
* @params
* @returns 
*/
function incompbeta(a, b, x){
     
    if(x == 0){
        return 0;
    }
    else if (x == 1){
        return 1;
    }else
    {
        let lbeta = lgamma(a+b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1-x)
        if (x < (a+1) / (a+b+2)){
            return Math.exp(lbeta) * contfractbeta(a, b, x) / a;
        }
        else{
            return 1 - Math.exp(lbeta) * contfractbeta(b, a, 1-x) / b;
        }
    }
}

/**
 * contfractbeta() evaluates the continued fraction form of the incomplete Beta function; incompbeta().  
 *   (Code translated from: Numerical Recipes in C.)
 * 
 * @param {*} a 
 * @param {*} b 
 * @param {*} x 
 * @param {*} ITMAX 
 * @returns 
 */
function contfractbeta(a,b,x, ITMAX = 1000){
     
    let EPS = 3.0e-7;
    let az = 1.0;
    let am = 1.0;
    let bm = 1.0;
    let qab = a + b
    let qap = a + 1.0
    let qam = a-1.0
    let bz = 1.0 - qab*x/qap
     
    //for i in range(ITMAX):
    for(let i =0; i<= ITMAX; i++){
        let em = parseFloat(i+1)
        let tem = em + em
        let d = em*(b-em)*x/((qam+tem)*(a+tem))
        let ap = az + d*am
        let bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
        let app = ap+d*az
        let bpp = bp+d*bz
        let aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.0
        if (Math.abs(az-aold)<(EPS * Math.abs(az))){
            return az
        }
    }
    return az // Return the last computed value if max iterations reached
}

/**
 * Evaluates factorial of a number
 * 
 * @param {Number} xf - Integer number 
 * @returns factorial of the number
 */
function factorial(xf) {
    if ((xf == 0) || (xf == 1)) return 1;
    else {
        let result = (xf * factorial(xf - 1))
        return result
    }
}

/**
 * Computes the natural logarithm of the Gamma function using Lanczos Approximation.
 * This avoids overflow for large values.
 *
 * @param {Number} x - Input value
 * @returns {Number} - Log Gamma function value
 */
function lgamma(x) {
    if (x < 0) return NaN; // Gamma is undefined for negative integers

    const g = 7;
    const coefficients = [
        0.99999999999980993, 
        676.5203681218851, 
        -1259.1392167224028, 
        771.32342877765313, 
        -176.61502916214059, 
        12.507343278686905, 
        -0.13857109526572012, 
        9.9843695780195716e-6, 
        1.5056327351493116e-7
    ];

    let z = x - 1;
    let x_sum = coefficients[0];
    
    for (let i = 1; i < coefficients.length; i++) {
        x_sum += coefficients[i] / (z + i);
    }

    let t = z + g + 0.5;
    return Math.log(Math.sqrt(2 * Math.PI)) + (z + 0.5) * Math.log(t) - t + Math.log(x_sum);
}


function t_test_1_sample(mean, m, s, n) {
    if (s == 0) s = 1;
    var t = ((mean - m) / s) * Math.sqrt(n)
    var p = 1.0 - TdistributionCDF(Math.abs(t), (n - 1))
    return p
}

function t_test_2_samples(m1, s1, n1, m2, s2, n2) {
    if (s1 == 0) s1 = 1;
    if (s2 == 0) s2 = 1;
    var t = (m1 - m2) / Math.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2);
    var df = ((s1 ** 2 / n1 + s2 ** 2 / n2) ** 2 * (n1 - 1) * (n2 - 1)) /
        ((s1 ** 4 * (n2 - 1)) / n1 ** 2 + (s2 ** 4 * (n1 - 1)) / n2 ** 2);

    var p = 1.0 - TdistributionCDF(Math.abs(t), parseInt(df + 0.5))

    return p
}

export default {TdistributionCDF, lgamma, t_test_1_sample, t_test_2_samples};