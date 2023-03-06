function solveTridiagonal (n, a, b, c, x) {
  var i, fac;

  // Eliminate:
  for (i = 1; i < n; i++) {
    if (b[i - 1] === 0) {
      console.log('tridiagonalSolve: failed due to lack of diagonal dominance');
      return false;
    }
    fac = a[i] / b[i - 1];
    b[i] -= fac * c[i - 1];
    x[i] -= fac * x[i - 1];
  }

  // Back-substitute:
  if (b[n - 1] === 0) {
    console.log('tridiagonalSolve: failed due to singular matrix');
    return false
  }
  x[n - 1] /= b[n - 1];
  for (i = n - 2; i >= 0; i--) {
    if (b[i] === 0) {
      console.log('tridiagonalSolve: failed due to singular matrix');
      return false;
    }
    x[i] = (x[i] - c[i] * x[i + 1]) / b[i];
  }

  return true;
}

function getScalarFieldID(r, sf,
  psi)
{
  NR = r.length;

  // radial differences
  var dr = Array(NR);
  dr[0] = r[1] - r[0];
  dr[NR-1] = r[NR-1] - r[NR-2];
  for(var i=1; i<NR-1; i++)
    dr[i] = (r[i+1] - r[i-1])/2.0;

  // radial derivative of scalar field profile
  var Phi = Array(NR).fill(0); // assume zero derivative at boundaries
  for(var i=1; i<NR-1; i++)
    Phi[i] = (sf[i+1] - sf[i-1]) / 2 / dr[i];

  // Set the main diagonal
  var main_diag = Array(NR);
  for(i=0; i<NR; i++)
    main_diag[i] = Math.PI*dr[i]*dr[i]*Phi[i]*Phi[i] - 2.0;

  // Update the first element of the main diagonal
  main_diag[0] += 1 - dr[0]/r[0]

  // Update the last element of the main diagonal
  main_diag[NR-1] += - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

  // Set the upper diagonal, ignoring the last point in the r array
  var upper_diag = Array(NR).fill(0);
  for(i=0; i<NR-1; i++)
    upper_diag[i] = 1 + dr[i]/r[i];

  // Set the lower diagonal, start counting the r array at the second element
  var lower_diag = Array(NR).fill(0);
  for(i=1; i<NR; i++)
    lower_diag[i] = 1 - dr[i]/r[i];

  // Change the last term in the lower diagonal to its correct value
  lower_diag[NR-1] = 2

  s = Array(NR).fill(0);
  // Update the last entry of the vector s
  s[NR-1] = - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

  // this replaces the RHS with the solution (psi)
  solveTridiagonal(NR, lower_diag, main_diag, upper_diag, s)

  // set the metric fields
  for(i=0; i<NR; i++)
    psi[i] = s[i];
  //console.log("Returning Phi Hopefully!", n) 
  //return Phi;
}

var NR = 30000;
var r = Array(NR),
    sf = Array(NR);
var psi = Array(NR),
    psi4 = Array(NR),
    alpha = Array(NR);

//Read in r0, phi0, and sigma as arguments (they define the gaussian)
function makeInitials(phi0=0.4, r0=0.0, sigma=1){
// coordinates
// below is for uniform radial spacing,
// check nrpy code for e.g. sinh spherical. More complicated...
// maybe ok here to just use more points.
    console.log('makeInitials:',phi0)
    console.log('makeInitials:',r0)
    console.log('makeInitials:',sigma)

    var NR = 30000,
        rmax = 12.0;
    for(var i=0; i<NR; i++)
      r[i] = (i+0.5)/NR*rmax;

// scalar field profile
    for(var i=0; i<NR; i++)
      sf[i] = phi0*Math.exp( -(r[i]-r0)*(r[i]-r0)/sigma/sigma );

    getScalarFieldID(r, sf, psi);
    
    console.log('Generating Psi4 and Alpha! Test!');
    for(var i=0; i<NR; i++)
      psi4[i] = Math.pow(psi[i],4);
    for(var i=0; i<NR; i++)
        alpha[i] = Math.pow(psi[i],-2);
    
    
}