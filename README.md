# Numerické metody

## Derivace

### Dopředná

```{r}
ForwardDiff <- function(f,x,h=1e-6) {
	return((f(x+h)-f(x))/h)
}
```

### Centrální diference

```{r}
CentralDiff <- function(f,x,h=1e-6) {
	return((f(x+h)-f(x-h))/(2*h))
}
```

## Integrály

### Pravidlo středního bodu

```{r}
MidpointRule <- function(f,a,b,n=1e6) {
	h <- (b-a)/n
	return(h*sum(f(a+h*(1:n)-h/2)))
}
```

### Simpsonovo pravidlo

```{r}
SimpsonRule <- function(f, a, b, n=1){
  h <- (b-a)/n
  hpul <- h/2
  suma <- f(a) + f(b)
  xs <- h*(1:n) + a - hpul
  suma <- suma + 4*sum(f(xs))
  if(n > 1){
    xl <- (xs + hpul)[-n]
    suma <- suma + 2*sum(f(xl))
  }
  return(suma*h/6)
}
```

### Monte Carlo

```{r}
GeomMethod <- function(f,a,b,h,n=1e6){
	x <- runif(n,a,b)
	y <- runif(n,0,h)
	return(h*(b-a)*sum(f(x) > y)/n)
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
Lagrange <- function(xa, x, y){
  n <- length(x)
  suma <- 0
  for(i in 1:n){
    nasobic <- 1
    for(j in 1:n){
      if(j != i) nasobic <- nasobic*(xa-x[j])/(x[i]-x[j])
    }
    suma <- suma + y[i]*nasobic
  }
  return(suma)
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

<!-- Spočte koeficienty aproximačního polynomu stupně $n - 1$ pro zadané body $x$, $y$. -->

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

### Bisekce

```{r}
BisecRoot <- function(f,a,b){
	fa <- f(a)
	fb <- f(b)
	if(fa*fb<0) {
		repeat {
			c <- (a+b)/2
			if(c==a | c==b) return(c)
			fc <- f(c)
			if(fa*fc<0) {
				b <- c
				fb <- fc
			} else {
				a <- c
				fa <- fc
			}
		}
	} else {
		stop("f(a)*f(b) < 0 not satisfied")
	}
}
```

<!-- Nalezne kořen spojité funkce $f(x)$ na intervalu $(a,b)$. $f(a)$ a $f(b)$ se musí lišit znaménkem, jinak není zaručeno, že se v daném intervalu kořen nachází. -->

### Newtonova metoda

```{r}
NewtonRoot <- function(f,x0,tol=1e-6) {
	x <- x0
    h <- tol
	repeat{
        dx <- 2*h*f(x)/(f(x+h)-f(x-h))
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
