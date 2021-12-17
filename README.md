# ndfwtree

N-dimensional Fenwick trees supporting different combinations of point/range updates/queries with dimensions set at compile-time. Requires fancy C++17 stuff like [fold expressions](https://en.cppreference.com/w/cpp/language/fold).

Not thoroughly tested. May be usable for competitive programming. Idk.

## Common

| Template Parameter | Description                                   | Default              |
| ------------------ | --------------------------------------------- | -------------------- |
| `typename T`       | type queried/updated, usually `int`           |                      |
| `T Id`             | identity value                                | `0`                  |
| `typename F`       | binary function                               | `std::plus<T>`       |
| `typename G`       | inverse binary function                       | `std::minus<T>`      |
| `typename H`       | binary function (repeated application of `F`) | `std::multiplies<T>` |
| `std::size_t N`    | number of dimensions                          | `1`                  |

`szarraytype` is a member type alias of `std::array<std::size_t, dimen>` where `dimen` is the number of dimensions. Size constructors and most member functions have an overload that takes a `szarraytype` instead of `dimen` number of `std::size_t`s.

Space complexity assumes `F`, `G`, and `H` are constant time (although they might not be). It also assumes the compiler will inline the recursive template instantiations (GCC seems like it does) and so does not count stack space.

## Range-Query + Point-Update

```cpp
template<typename T, T Id = 0, typename F = std::plus<T>, typename G = std::minus<T>, std::size_t N = 1>
class ndfwtree;
```

**Complexity**

- Space: O(d+n^d)

### Constructors

```cpp
template<typename F1 = F, typename G1 = G>
explicit ndfwtree(std::size_t sz_0, ..., std::size_t sz_n, F1&& f = F(), G1&& g = G());
template<typename F1 = F, typename G1 = G>
explicit ndfwtree(const szarraytype& szs, F1&& f = F(), G1&& g = G());
```

Size constructors. Constructs Fenwick tree from the size of each dimension, initialized to `Id`.

```cpp
template<typename V, typename F1 = F, typename G1 = G>
ndfwtree(V&& v, F1&& f = F(), G1&& g = G());
template<typename V, typename F1 = F, typename G1 = G>
ndfwtree(V&& v, std::size_t sz_0, ..., std::size_t sz_n, F1&& f = F(), G1&& g = G());
```

Flat vector build constructors. Builds Fenwick tree using `v`. Both constructors only participate in overload resolution if `v` is a reference to (cv-qualified) `std::vector<T>`. The first constructor requires dimensions to be 1 and the second requires dimensions to be >=2.

```cpp
template<typename F1 = F, typename G1 = G>
ndfwtree(const vector<...<vector<T>>...>& v, F1&& f = F(), G1&& g = G())
```

N-dimensional vector build constructor. Builds Fenwick tree using `v`. Requires dimensions to be >= 2.

**Build Complexity**

- Time: O(d*n^d)
- Aux space: O(d)

### Member Functions

```cpp
T prefixqry(std::size_t r_0, ..., std::size_t r_n) const;
T prefixqry(const szarraytype& rs) const;
```

Query `[0, r)` in each dimension.

**Complexity**

- Time: O(log^d(n))

- Aux space: O(d)

```cpp
T rangeqry(std::size_t l_0, ..., std::size_t l_n, std::size_t r_0, ..., std::size_t r_n) const;
T rangeqry(const szarraytype& ls, const szarraytype& rs) const;
```

Query `[l, r)` in each dimension.

**Complexity**

- Time: O((2*log(n))^d)

- Aux space: O(d)

```cpp
void pointupd(std::size_t i_0, ..., std::size_t i_n, T delta);
void pointupd(const szarraytype& is, T delta);
void pointupdinv(std::size_t i_0, ..., std::size_t i_n, T delta);
void pointupdinv(const szarraytype& is, T delta);
```

Update `i` in each dimension.

**Complexity**

- Time: O(log^d(n))

- Aux space: O(d)

```cpp
template<size_t I = 0>
size_t size() const;
F f() const;
G g() const;
```

## Point-Query + Range-Update

```cpp
template<typename T, T Id = 0, typename F = plus<T>, typename G = minus<T>, size_t N = 1>
class ndfwtreepr;
```

**Complexity**

- Space: O(d+n^d)

### Constructors

```cpp
template<typename F1 = F, typename G1 = G>
explicit ndfwtreepr(std::size_t sz_0, ..., std::size_t sz_n, F1&& f = F(), G1&& g = G());
template<typename F1 = F, typename G1 = G>
explicit ndfwtreepr(const szarraytype& szs, F1&& f = F(), G1&& g = G());
```

Size constructors. Constructs Fenwick tree from the size of each dimension, initialized to `Id`.

```cpp
template<typename V, typename F1 = F, typename G1 = G>
ndfwtreepr(V&& v, F1&& f = F(), G1&& g = G());
template<typename V, typename F1 = F, typename G1 = G>
ndfwtreepr(V&& v, std::size_t sz_0, ..., std::size_t sz_n, F1&& f = F(), G1&& g = G());
```

Flat vector build constructors. Builds Fenwick tree using `v`. Both constructors only participate in overload resolution if `v` is a reference to (cv-qualified) `std::vector<T>`. The first constructor requires dimensions to be 1 and the second requires dimensions to be >=2.

```cpp
template<typename F1 = F, typename G1 = G>
ndfwtreepr(const vector<...<vector<T>>...>& v, F1&& f = F(), G1&& g = G());
```

N-dimensional vector build constructor. Builds Fenwick tree using `v`. Requires dimensions to be >= 2.

**Build Complexity**

- Time: O(d\*n\^d+(2*n)^d)
- Aux space: O(d)

### Member Functions

```cpp
T pointqry(std::size_t i_0, ..., std::size_t i_n) const;
T pointqry(const szarraytype& is) const;
```

Query `i` in each dimension.

**Complexity**

- Time: O(log^d(n))
- Aux space: O(d)

```cpp
void suffixupd(std::size_t l_0, ..., std::size_t l_n, T delta);
void suffixupd(const szarraytype& ls, T delta);
void suffixupdinv(std::size_t l_0, ..., std::size_t l_n, T delta);
void suffixupdinv(const szarraytype& ls, T delta);
```

Update `[l, sz)` in each dimension.

**Complexity**

- Time: O(log^d(n))
- Aux space: O(d)

```cpp
void rangeupd(std::size_t l_0, ..., std::size_t ls_n, std::size_t r_0, ..., std::size_t r_n, T delta);
void rangeupd(const szarraytype& ls, const szarraytype& rs, T delta);
void rangeupdinv(std::size_t l_0, ..., std::size_t l_n, std::size_t r_0, ..., std::size_t r_n, T delta) ;
void rangeupdinv(const szarraytype& ls, const szarraytype& rs, T delta);
```

Update `[l, r)` in each dimension.

**Complexity**

- Time: O((2*log(n))^d)
- Aux space: O(d)

```cpp
template<size_t I = 0>
size_t size() const;
F f() const;
G g() const;
```

## Range-Query + Range-Update

```cpp
template<typename T, T Id = 0, typename F = plus<T>, typename G = minus<T>, typename H = multiplies<T>, size_t N = 1>
class ndfwtreerr;
```

**Complexity**

- Space: O(d+(2*n)^d)

### Constructors

```cpp
template<typename F1 = F, typename G1 = G, typename H1 = H>
explicit ndfwtreerr(std::size_t sz_0, ..., std::size_t sz_n, F1&& f = F(), G1&& g = G(), H1&& h = H());
template<typename F1 = F, typename G1 = G, typename H1 = H>
explicit ndfwtreerr(const szarraytype& szs, F1&& f = F(), G1&& g = G(), H1&& h = H());
```

Size constructors. Constructs Fenwick tree from the size of each dimension, initialized to `Id`.

```cpp
template<typename F1 = F, typename G1 = G, typename H1 = H>
explicit ndfwtreerrimpl(const vector<...<vector<T>>...>& v, F1&& f = F(), G1&& g = G(), H1&& h = H());
```

N-dimensional vector build constructor. Builds Fenwick tree using `v`.

**Build Complexity**

- Time: O(d\*(2\*n)\^d+(4*n)^d)
- Aux space: O(d)

### Member Functions

```cpp
T prefixqry(std::size_t r_0, ..., std::size_t r_n) const;
T prefixqry(const szarraytype& rs) const;
```

Query `[0, r)` in each dimension.

**Complexity**

- Time: O((2*log(n))^d)
- Aux space: O(d)

```cpp
T rangeqry(std::size_t l_0, ..., std::size_t l_n, std::size_t r_0, ..., std::size_t r_n) const;
T rangeqry(const szarraytype& ls, const szarraytype& rs) const;
```

Query `[l, r)` in each dimension.

**Complexity**

- Time: O((4*log(n))^d)
- Aux space: O(d)

```cpp
void suffixupd(std::size_t l_0, ..., std::size_t l_n, T delta);
void suffixupd(const szarraytype& ls, T delta);
void suffixupdinv(std::size_t l_0, ..., std::size_t l_n, T delta);
void suffixupdinv(const szarraytype& ls, T delta);
```

Update `[l, sz)` in each dimension.

**Complexity**

- Time: O((2*log(n))^d)
- Aux space: O(d)

```cpp
void rangeupd(std::size_t l_0, ..., std::size_t ls_n, std::size_t r_0, ..., std::size_t r_n, T delta);
void rangeupd(const szarraytype& ls, const szarraytype& rs, T delta);
void rangeupdinv(std::size_t l_0, ..., std::size_t ls_n, std::size_t r_0, ..., std::size_t r_n, T delta);
void rangeupdinv(const szarraytype& ls, const szarraytype& rs, T delta);
```

Update `[l, r)` in each dimension.

**Complexity**

- Time: O((4*log(n))^d)
- Aux space: O(d)

```cpp
template<size_t I = 0>
size_t size() const;
F f() const;
G g() const;
H h() const;
```

## Aliases

```cpp
template<typename T>
using fwtree = ndfwtree<T>;
template<typename T>
using fwtreepr = ndfwtreepr<T>;
template<typename T>
using fwtreerr = ndfwtreerr<T>;
template<typename T>
using fwtreemul = ndfwtree<T, 1, multiplies<T>, divides<T>>;
template<typename T>
using fwtreeprmul = ndfwtreepr<T, 1, multiplies<T>, divides<T>>;
template<typename T>
using fwtreerrmul = ndfwtreerr<T, 1, multiplies<T>, divides<T>, powfn<T>>;
template<typename T>
using fwtree2d = ndfwtree<T, 0, plus<T>, minus<T>, 2>;
template<typename T>
using fwtreepr2d = ndfwtreepr<T, 0, plus<T>, minus<T>, 2>;
template<typename T>
using fwtreerr2d = ndfwtreerr<T, 0, plus<T>, minus<T>, multiplies<T>, 2>;
template<typename T>
using fwtreemul2d = ndfwtree<T, 1, multiplies<T>, divides<T>, 2>;
template<typename T>
using fwtreeprmul2d = ndfwtreepr<T, 1, multiplies<T>, divides<T>, 2>;
template<typename T>
using fwtreerrmul2d = ndfwtreerr<T, 1, multiplies<T>, divides<T>, powfn<T>, 2>;
```

`powfn` is a functor for the [power](https://en.cppreference.com/w/cpp/numeric/math/pow) operation.

# Resources

1. https://cp-algorithms.com/data_structures/fenwick.html
2. https://arxiv.org/pdf/1311.6093v4.pdf (range-range)
3. https://acm.timus.ru/problem.aspx?space=1&num=1470 (practice problem)
