# Přehled vybraných numerických metod

## Derivace

### Dopředná

```{r}
d <- function(f,x,h) {
	return((f(x+h)-f(x))/h)
}
```

### Centrální diference

```{r}
d <- function(f,x,h) {
	return((f(x+h)-f(x-h))/(2*h))
}
```

## Integrály

### Midpoint rule

```{r}
MidpointRule <- function(f,a,b,n) {
	h <- (b-a)/n
	return(h*sum(f(a+h*(1:n)-h/2)))
}
```

### Monte Carlo

```{r}
GeomMethod <- function(f,a,b,h,n){
	x <- runif(n,a,b)
	y <- runif(n,0,h)
	return(h*(b-a)*sum(f(x) > y)/n)
}
```

```{r}
AveMethod <- function(f,a,b,n){
	return((b-a)*mean(f(runif(n,a,b))))
}
```

## Soustavy rovnic

### Gaussova eliminace

```{r}
GaussPivot <- function(A,b) {
	Ab <- cbind(A, b)
	n <- length(b)
	for(k in 1:(n-1)) {
		pivot <- which.max(abs(Ab[k:n,k])) + k - 1
		if(pivot != k){
			j <- k:(n+1)
			pom <- Ab[k,j]
			Ab[k,j] <- Ab[pivot,j]
			Ab[pivot,j] <- pom
		}
		for (i in (k+1):n) {
			j <- (k+1):(n+1)
			nasobek <- Ab[i,k] / Ab[k,k]
			Ab[i,j] <- Ab[i,j] - nasobek*Ab[k,j]
		}
	}
	x <- b
	x[n] <- Ab[n,n+1] / Ab[n,n]
	for(i in (n-1):1){
		j <- (i+1):n
		x[i] <- (Ab[i,n+1]-sum(Ab[i,j]*x[j]))/Ab[i,i]
	}
	return(x)
}
```

## Interpolace

<!--
## Newton
-->

### Lagrange

```{r}
l <- function(xa, x, j){
  res <- 1
  for(i in 1:length(x)){
    if(i != j) res <- res*(xa-x[i])/(x[j]-x[i])
  }
  return(res)
}
Lagrange <- function(xa, x, y){
  res <- 0
  for(i in 1:length(x)) res <- res+y[i]*l(xa, x, i)
  return(res)
}
```

## Aproximace

### LSA

```{r}
LSA <- function(x, y, n){
  X <- matrix(1, nrow=n, ncol=length(x))
  for(i in 2:n) X[i, ] <- X[i-1, ]*x
  return(c(solve(X%*%t(X), X%*%y)))
}
```

Spočte koeficienty aproximačního polynomu stupně $n - 1$ pro zadané body $x$, $y$.

## Diferenciální rovnice

### Eulerova metoda

```{r}
EulerStep <- function(f,x,y,h) {
	return(y+h*f(x,y))
}
```

### RK4

```{r}
RK4 <- function(f,x,y,h){
  hhalf <- 0.5*h
  k1 <- f(x,y)
  k2 <- f(x+hhalf, y+hhalf*k1)
  k3 <- f(x+hhalf, y+hhalf*k2)
  k4 <- f(x+h, y+h*k3)
  return(y+h*(k1+2*(k2+k3)+k4)/6)
}
```

### Metoda prosté iterace

```{r}
FixedPointStep <- function(f, h, x0, y0){
  hpul <- 0.5*h
  x1 <- x0+h
  y1 <- y0+h*f(x0,y0)
  for(i in 1:10) y1 <- y0+hpul*(f(x1,y1)+f(x0,y0))
  return(y1)
}
```

## Hledání kořenů

### Newtonova metoda

```{r}
NewtonRoot <- function(f,fd,x0,tol=1e-6) {
	x <- x0
	repeat{
		dx <- f(x)/fd(x)
		if(abs(dx) < tol) return(x)
		x <- x - dx
	}
}
```

## Polynomy

### Horner

```{r}
Horner <- function(a,x){
	n <- length(a)
	y <- a[n]
	for(i in (n-1):1) y <- y*x+a[i]
	return(y)
}
```

### Chebyshev

```{r}
ChebyCoef <- function(n){
	a0 <- numeric(n)
	a0[1] <- 1
	if(n==1) return(a0)
	a1 <- numeric(n)
	a1[2] <- 1
	if(n==2) return(a1)
	for(i in 3:n){
		a <- 2*c(0, a1[-n])-a0
		a0 <- a1
		a1 <- a
	}
	return(a)
}
```
