Omegas = function(bspline.basis, ss)   {
    require(fda)
    k.n = length(bspline.basis$params)
    k = c(rep(bspline.basis$rangeval[1],4), bspline.basis$params, rep(bspline.basis$rangeval[2],4))
    hh = rep(NA, k.n+7)
    for (i in 1:(k.n+7))    hh[i] = (k[i+1]-k[i])/(6-2*ss)

    # Table 1 of Wand and Ormerod (2008)
    w.mat = rbind(c(1, NA, NA, NA, NA, NA, NA),
                  c(1/3, 4/3, 1/3, NA, NA, NA, NA),
                  c(14/45, 64/45, 8/15, 64/45, 14/45, NA, NA),
                  c(41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140))
    x = w = rep(NA, (7-2*ss)*(k.n+7))
    for (l in 1:(k.n+7))
        for (ll in 0:(6-2*ss))   {
            x[(7-2*ss)*(l-1)+ll+1]= k[l]+ll*hh[l]
            w[(7-2*ss)*(l-1)+ll+1]= hh[l]*w.mat[4-ss,ll+1]
        }
    B = eval.basis(x, bspline.basis, ss)
    t(B) %*% diag(w) %*% B
}

