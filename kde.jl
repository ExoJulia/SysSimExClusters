using Distributions
using QuadGK

function kde_kth_neighbor_1d(xout::Vector{T}, xin::Vector{T}; k::Integer = 10, sigma::Real = 2.0) where T <: Real
   issorted(xin) ? kde_kth_neighbor_1d_sorted(xout,xin,k=k,sigma=sigma)  :
                   kde_kth_neighbor_1d_unsorted(xout,xin,k=k,sigma=sigma)  
end

function kde_kth_neighbor_1d_unsorted(xout::Vector{T}, xin::Vector{T}; k::Integer = 10, sigma::Real = 2) where T <: Real
   @assert(length(xin)>=k+1)
   @assert issorted(xout)
   yout = similar(xout)
   ilo = 1
   ihi = k
   perm = sortperm(xin)
   #print("# ilo,ihi = ",ilo,", ", ihi, " length(xin)=", length(xin), " len(xout)=",length(xout) )
   #println( " perm[ilo,ihi,ihi+1] = ", perm[ilo],", ", perm[ihi], ", ", perm[ihi+1])
   dlo = abs(xout[1]-xin[perm[ilo]])
   dhi = abs(xin[perm[ihi]]-xout[1])
   dstar = abs(xin[perm[ihi+1]]-xout[1])
   for (i,x) in enumerate(xout)
       while dstar<dlo && ihi<length(xin) 
          ilo += 1 
          ihi += 1
          dlo = abs(x-xin[perm[ilo]])
          dhi = dstar
          if ihi+1<=length(perm)
            dstar = abs(xin[perm[ihi+1]]-x)
          end
       end 
       #width = (xin[perm[ihi]]-xin[perm[ilo]])*sigma/k
       width = max(dlo,dhi)*sigma/k
       distrib = Distributions.Normal(x,width)
       yout[i] = sum(pdf(distrib,xin[perm[ilo:ihi]]))/k
   end
   return yout
end

function kde_kth_neighbor_1d_sorted(xout::Vector{T}, xin::Vector{T}; k::Integer = 10, sigma::Real = 2) where T <: Real
   @assert(length(xin)>=k+1)
   #@assert issorted(xin)
   yout = similar(xout)
   ilo = 1
   ihi = k
   dlo = abs(xout[1]-xin[ilo])
   dhi = abs(xin[ihi]-xout[1])
   dstar = abs(xin[ihi+1]-xout[1])
   for (i,x) in enumerate(xout)
       while dstar<dlo && ihi<length(xin) 
          ilo += 1 
          ihi += 1
          dlo = abs(x-xin[ilo])
          dhi = dstar
          if ihi+1<=length(perm)
             dstar = abs(xin[ihi+1]-x)
          end
       end 
       #width = (xin[ihi]-xin[ilo])*sigma/k
       width = max(dlo,dhi)*sigma/k
       distrib = Distributions.Normal(x,width)
       yout[i] = sum(pdf(distrib,xin[ilo:ihi]))/k
   end
   return yout
end

function kl_integral_arg(q::T, p::T; epsilon::T = 1e-10) where T <: Real
  if q > epsilon
     return p*log(p/q)
  else
     return p
  end
end

function hellinger_integral_arg(p::T, q::T) where T <: Real
  sqrt(p*q)
end


function calc_kl_distance_ab(x1::Vector{T}, x2::Vector{T}, a::T, b::T; n::Integer = 100, k::Integer = max(min((length(x1) - 1), (length(x2) - 1), 20), 1)) where T <: Real
  @assert b>a
  @assert n>=2
  if(length(x1)<3 || length(x2)<3) return 0.0  end
  xgrid = collect(linspace(a,b,n+1))
  y1 = kde_kth_neighbor_1d(xgrid,x1,k=k)
  y2 = kde_kth_neighbor_1d(xgrid,x2,k=k)
  integral = 0.5*(kl_integral_arg(y1[1],y2[1])+kl_integral_arg(y1[end],y2[end]))
  integral += sum([kl_integral_arg(y1[i],y2[i]) for i in 2:(length(y1)-1) ])
  integral /= n
end

function calc_kl_distance(x1::Vector{T}, x2::Vector{T}; n::Integer = 100, k::Integer = max(min((length(x1) - 1), (length(x2) - 1), 20), 1)) where T <: Real
  @assert b>a
  @assert n>=2
  n1 = length(x1)
  n2 = length(x2)
  xmin = minimum(x1,x2)
  xmax = maximum(x1,x2)
  padding = k*(xmax-xmin)/(n1+n2)
  calc_kl_distance_ab(x1,x2,xmin,xmax,n=n,k=k)
end



