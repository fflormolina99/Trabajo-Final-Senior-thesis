Module randomize

IMPLICIT none 

Contains

function mzran()
integer mzran,mzranset,i,j,k,n,is,js,ks,ns
save i,j,k,n
data i,j,k,n/521288629,362436069,16163801,1131199299/
mzran = i-k
if(mzran.lt.0) mzran = mzran + 2147483579
i = j
j = k
k = mzran
n = 69069*n + 1013904243
mzran = mzran + n
return
entry mzranset(is,js,ks,ns)
i = 1+iabs(is)
j = 1+iabs(js)
k = 1+iabs(ks)
n = ns
mzranset = n
return
end function mzran

function random()
real random
random = .5+.2328306E-9*mzran()
return
end function random

end module randomize
