osplinepen2d = function(bbasis1, bbasis2)   {
    omega20 = Omegas(bbasis1, 2) %x% Omegas(bbasis2, 0)
    omega11 = Omegas(bbasis1, 1) %x% Omegas(bbasis2, 1)
    omega02 = Omegas(bbasis1, 0) %x% Omegas(bbasis2, 2)
    omega20 + 2 * omega11 + omega02
}
