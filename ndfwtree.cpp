#include <bits/stdc++.h>
using namespace std;

// https://cp-algorithms.com/data_structures/fenwick.html

namespace detail {

// compiler gods please inline
template<typename T, size_t N>
struct ndvector;
template<typename T, size_t N>
using ndvectortype = typename ndvector<T, N>::type;
template<typename T, size_t N>
struct ndvector { using type = vector<ndvectortype<T, N - 1>>; };
template<typename T>
struct ndvector<T, 1> { using type = vector<T>; };
template<typename V>
struct ndvectordimen : integral_constant<size_t, 0> {};
template<typename V>
constexpr size_t ndvectordimenvalue = ndvectordimen<V>::value;
template<typename V>
struct ndvectordimen<vector<V>>
    : integral_constant<size_t, ndvectordimenvalue<V> + 1> {};
template<size_t I, typename V, size_t N>
inline void ndvectorsizeshelper(const V& v, array<size_t, N>& szs) {
    szs[I] = v.size();
    if constexpr (I < N - 1)
        ndvectorsizeshelper<I + 1>(v[0], szs);
}
template<typename V>
inline auto ndvectorsizes(const V& v) {
    array<size_t, ndvectordimenvalue<V>> szs;
    ndvectorsizeshelper<0>(v, szs);
    return szs;
}

template<size_t>
using alwaysszt = size_t;
template<typename T>
using removecvreftype = remove_cv_t<remove_reference_t<T>>;

template<typename T, typename U>
void copyndvector(vector<T>& bit, const U& u, size_t idx = 0) {
    if constexpr (is_same_v<U, T>) {
        bit[idx] = u;
    } else {
        idx *= u.size();
        for (size_t i = 0; i < u.size(); ++i)
            copyndvector(bit, u[i], idx + i);
    }
}

template<typename T>
inline T lsb(T x) {
    return x & -static_cast<uintmax_t>(x);
}

template<size_t I, typename T, size_t N, typename F>
void buildbithelper(
    vector<T>& bit, const array<size_t, N>& szs, const F& f,
    size_t lhs, size_t rhs)
{
    if constexpr (I == N) {
        bit[lhs] = f(bit[lhs], bit[rhs]);
    } else {
        for (size_t i = 0; i < szs[I]; ++i)
            buildbithelper<I + 1>(
                bit, szs, f, lhs * szs[I] + i, rhs * szs[I] + i);
    }
}
template<size_t I, typename T, size_t N, typename F>
void buildbit(
    vector<T>& bit, const array<size_t, N>& szs, const F& f, size_t idx = 0)
{
    idx *= szs[I];
    if constexpr (I < N - 1)
        for (size_t i = 0; i < szs[I]; ++i)
            buildbit<I + 1>(bit, szs, f, idx+ i);
    for (size_t i = 0; i < szs[I]; ++i) {
        size_t j = i + lsb(i + 1);
        if (j < szs[I])
            buildbithelper<I + 1>(bit, szs, f, idx + j, idx + i);
    }
}
template<typename T, size_t K, size_t N, typename F, size_t... Idxs>
inline void buildbitshelper1(
    vector<array<T, K>>& bits, const array<size_t, N>& szs, const F& f,
    size_t lhs, size_t rhs, index_sequence<Idxs...>)
{
    ((bits[lhs][Idxs] = f(bits[lhs][Idxs], bits[rhs][Idxs])), ...);
}
template<size_t I, typename T, size_t K, size_t N, typename F>
void buildbitshelper(
    vector<array<T, K>>& bits, const array<size_t, N>& szs, const F& f,
    size_t lhs, size_t rhs)
{
    if constexpr (I == N) {
        buildbitshelper1(bits, szs, f, lhs, rhs, make_index_sequence<K>());
    } else {
        lhs *= szs[I];
        rhs *= szs[I];
        for (size_t i = 0; i < szs[I]; ++i)
            buildbitshelper<I + 1>(bits, szs, f, lhs + i, rhs + i);
    }
}
template<size_t I, typename T, size_t K, size_t N, typename F>
void buildbits(
    vector<array<T, K>>& bits, const array<size_t, N>& szs, const F& f,
    size_t idx = 0)
{
    idx *= szs[I];
    if constexpr (I < N - 1)
        for (size_t i = 0; i < szs[I]; ++i)
            buildbits<I + 1>(bits, szs, f, idx + i);
    for (size_t i = 0; i < szs[I]; ++i) {
        size_t j = i + lsb(i + 1);
        if (j < szs[I])
            buildbitshelper<I + 1>(bits, szs, f, idx + j, idx + i);
    }
}

template<typename T, T Id, typename F, typename G, typename Seq>
class ndfwtreeimpl;
template<typename T, T Id, typename F, typename G, typename Seq>
class ndfwtreeprimpl;
template<typename T, T Id, typename F, typename G, typename H, typename Seq>
class ndfwtreerrimpl;

template<typename T, T Id, typename F, typename G, size_t... Idxs>
class ndfwtreeimpl<T, Id, F, G, index_sequence<Idxs...>> {
public:
    static constexpr size_t dimen = sizeof...(Idxs);
    static constexpr T identity = Id;
    using valuetype = T;
    using ftype = F;
    using gtype = G;
    using szarraytype = array<size_t, dimen>;
    template<typename F1 = F, typename G1 = G>
    explicit ndfwtreeimpl(alwaysszt<Idxs>... szs, F1&& f = F(), G1&& g = G())
        : szs_{szs...}, bit_((... * szs))
        , f_(forward<F1>(f)), g_(forward<G1>(g)) {}
    template<typename F1 = F, typename G1 = G>
    explicit ndfwtreeimpl(const szarraytype& szs, F1&& f = F(), G1&& g = G())
        : ndfwtreeimpl(
            szs, forward<F1>(f), forward<G1>(g),
            make_index_sequence<dimen>()) {}
    template<
        typename V, typename F1 = F, typename G1 = G,
        enable_if_t<
            dimen == 1 && is_same_v<removecvreftype<V>, vector<T>>>* = nullptr>
    ndfwtreeimpl(V&& v, F1&& f = F(), G1&& g = G())
        : szs_{v.size()}, bit_(forward<V>(v))
        , f_(forward<F1>(f)), g_(forward<G1>(g))
    {
        buildbit<0>(bit_, szs_, f_);
    }
    template<
        typename V, typename F1 = F, typename G1 = G,
        enable_if_t<
            dimen >= 2 && is_same_v<removecvreftype<V>, vector<T>>>* = nullptr>
    ndfwtreeimpl(V&& v, alwaysszt<Idxs>... szs, F1&& f = F(), G1&& g = G())
        : szs_{szs...}, bit_(forward<V>(v))
        , f_(forward<F1>(f)), g_(forward<G1>(g))
    {
        buildbit<0>(bit_, szs_, f_);
    }
    template<
        typename F1 = F, typename G1 = G, size_t N = dimen,
        enable_if_t<N >= 2>* = nullptr>
    ndfwtreeimpl(
        const ndvectortype<T, dimen>& v, F1&& f = F(), G1&& g = G())
        : ndfwtreeimpl(ndvectorsizes(v), forward<F1>(f), forward<G1>(g))
    {
        copyndvector(bit_, v);
        buildbit<0>(bit_, szs_, f_);
    }
    T prefixqry(alwaysszt<Idxs>... rs) const {
        return prefixqry(array{rs...});
    }
    T prefixqry(const szarraytype& rs) const {
        return prefixqryhelper<0>(rs);
    }
    T rangeqry(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs) const {
        return rangeqry(array{ls...}, array{rs...});
    }
    T rangeqry(const szarraytype& ls, const szarraytype& rs) const {
        return rangeqryhelper<0>(ls, rs);
    }
    void pointupd(alwaysszt<Idxs>... is, T delta) {
        pointupd(array{is...}, delta);
    }
    void pointupd(const szarraytype& is, T delta) {
        pointupdhelper<true, 0>(is, delta);
    }
    void pointupdinv(alwaysszt<Idxs>... is, T delta) {
        pointupdinv(array{is...}, delta);
    }
    void pointupdinv(const szarraytype& is, T delta) {
        pointupdhelper<false, 0>(is, delta);
    }
    template<size_t I = 0>
    size_t size() const { return szs_[I]; }
    F f() const { return f_; }
    G g() const { return g_; }
private:
    template<typename F1, typename G1, size_t... Idxs1>
    ndfwtreeimpl(
        const szarraytype& szs, F1&& f, G1&& g, index_sequence<Idxs1...>)
        : ndfwtreeimpl(szs[Idxs1]..., forward<F1>(f), forward<G1>(g)) {}
    template<size_t I>
    T prefixqryhelper(const szarraytype& rs, size_t idx = 0) const {
        if constexpr (I == dimen) {
            return bit_[idx];
        } else {
            idx *= szs_[I];
            T res = Id;
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                res = f_(res, prefixqryhelper<I + 1>(rs, idx + r - 1));
            return res;
        }
    }
    template<size_t I>
    T rangeqryhelper(
        const szarraytype& ls, const szarraytype& rs, size_t idx = 0) const
    {
        if constexpr (I == dimen) {
            return bit_[idx];
        } else {
            idx *= szs_[I];
            T res = Id;
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                res = f_(res, rangeqryhelper<I + 1>(ls, rs, idx + r - 1));
            for (size_t l = ls[I]; l > 0; l -= lsb(l))
                res = g_(res, rangeqryhelper<I + 1>(ls, rs, idx + l - 1));
            return res;
        }
    }
    template<bool B, size_t I>
    void pointupdhelper(const szarraytype& is, T delta, size_t idx = 0)
    {
        if constexpr (I == dimen) {
            if constexpr (B)
                bit_[idx] = f_(bit_[idx], delta);
            else
                bit_[idx] = g_(bit_[idx], delta);
        } else {
            idx *= szs_[I];
            for (size_t i = is[I]; i < szs_[I]; i += lsb(i + 1))
                pointupdhelper<B, I + 1>(is, delta, idx + i);
        }
    }
    szarraytype szs_;
    vector<T> bit_;
    F f_; G g_;
};

template<typename T, T Id, typename F, typename G, size_t... Idxs>
class ndfwtreeprimpl<T, Id, F, G, index_sequence<Idxs...>> {
public:
    static constexpr size_t dimen = sizeof...(Idxs);
    static constexpr T identity = Id;
    using valuetype = T;
    using ftype = F;
    using gtype = G;
    using szarraytype = array<size_t, dimen>;
    template<typename F1 = F, typename G1 = G>
    explicit ndfwtreeprimpl(alwaysszt<Idxs>... szs, F1&& f = F(), G1&& g = G())
        : szs_{szs...}, bit_((... * szs))
        , f_(forward<F1>(f)), g_(forward<G1>(g)) {}
    template<typename F1 = F, typename G1 = G>
    explicit ndfwtreeprimpl(const szarraytype& szs, F1&& f = F(), G1&& g = G())
        : ndfwtreeprimpl(
            szs, forward<F1>(f), forward<G1>(g),
            make_index_sequence<dimen>()) {}
    template<
        typename V, typename F1 = F, typename G1 = G,
        enable_if_t<
            dimen == 1 && is_same_v<removecvreftype<V>, vector<T>>>* = nullptr>
    ndfwtreeprimpl(V&& v, F1&& f = F(), G1&& g = G())
        : szs_{v.size()}, bit_(forward<V>(v))
        , f_(forward<F1>(f)), g_(forward<G1>(g))
    {
        szarraytype is;
        preprocess<0>(is);
        buildbit<0>(bit_, szs_, f_);
    }
    template<
        typename V, typename F1 = F, typename G1 = G,
        enable_if_t<
            dimen >= 2 && is_same_v<removecvreftype<V>, vector<T>>>* = nullptr>
    ndfwtreeprimpl(V&& v, alwaysszt<Idxs>... szs, F1&& f = F(), G1&& g = G())
        : szs_{szs...}, bit_(forward<V>(v))
        , f_(forward<F1>(f)), g_(forward<G1>(g))
    {
        szarraytype is;
        preprocess<0>(is);
        buildbit<0>(bit_, szs_, f_);
    }
    template<
        typename F1 = F, typename G1 = G, size_t N = dimen,
        enable_if_t<N >= 2>* = nullptr>
    ndfwtreeprimpl(
        const ndvectortype<T, dimen>& v, F1&& f = F(), G1&& g = G())
        : ndfwtreeprimpl(ndvectorsizes(v), forward<F1>(f), forward<G1>(g))
    {
        copyndvector(bit_, v);
        szarraytype is;
        preprocess<0>(is);
        buildbit<0>(bit_, szs_, f_);
    }
    T pointqry(alwaysszt<Idxs>... is) const {
        return pointqry(array{is...});
    }
    T pointqry(const szarraytype& is) const {
        return pointqryhelper<0>(is);
    }
    void suffixupd(alwaysszt<Idxs>... ls, T delta) {
        suffixupd(array{ls...}, delta);
    }
    void suffixupd(const szarraytype& ls, T delta) {
        suffixupdhelper<true, 0>(ls, delta);
    }
    void suffixupdinv(alwaysszt<Idxs>... ls, T delta) {
        suffixupdinv(array{ls...}, delta);
    }
    void suffixupdinv(const szarraytype& ls, T delta) {
        suffixupdhelper<false, 0>(ls, delta);
    }
    void rangeupd(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs, T delta) {
        rangeupd(array{ls...}, array{rs...}, delta);
    }
    void rangeupd(const szarraytype& ls, const szarraytype& rs, T delta) {
        rangeupdhelper<true, 0>(ls, rs, delta);
    }
    void rangeupdinv(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs, T delta) {
        rangeupdinv(array{ls...}, array{rs...}, delta);
    }
    void rangeupdinv(const szarraytype& ls, const szarraytype& rs, T delta) {
        rangeupdhelper<false, 0>(ls, rs, delta);
    }
    template<size_t I = 0>
    size_t size() const { return szs_[I]; }
    F f() const { return f_; }
    G g() const { return g_; }
private:
    template<typename F1, typename G1, size_t... Idxs1>
    ndfwtreeprimpl(
        const szarraytype& szs, F1&& f, G1&& g, index_sequence<Idxs1...>)
        : ndfwtreeprimpl(szs[Idxs1]..., forward<F1>(f), forward<G1>(g)) {}
    template<size_t I>
    void preprocess(szarraytype& is, size_t idx = 0) {
        if constexpr (I == dimen) {
            preprocesshelper<true, true, 0>(is, bit_[idx]);
        } else {
            idx *= szs_[I];
            for (is[I] = szs_[I]; is[I]-- > 0;)
                preprocess<I + 1>(is, idx + is[I]);
        }
    }
    template<bool B, bool Same, size_t I>
    void preprocesshelper(szarraytype& is, T delta, size_t idx = 0) {
        if constexpr (I == dimen) {
            // Same is true when idx is current element
            if constexpr (!Same) {
                if constexpr (B)
                    bit_[idx] = f_(bit_[idx], delta);
                else
                    bit_[idx] = g_(bit_[idx], delta);
            }
        } else {
            idx *= szs_[I];
            preprocesshelper<B, Same, I + 1>(is, delta, idx + is[I]);
            if (is[I] + 1 < szs_[I])
                preprocesshelper<!B, false, I + 1>(is, delta, idx + is[I] + 1);
        }
    }
    template<size_t I>
    T pointqryhelper(const szarraytype& is, size_t idx = 0) const {
        if constexpr (I == dimen) {
            return bit_[idx];
        } else {
            idx *= szs_[I];
            T res = Id;
            for (size_t i = is[I] + 1; i > 0; i -= lsb(i))
                res = f_(res, pointqryhelper<I + 1>(is, idx + i - 1));
            return res;
        }
    }
    template<bool B, size_t I>
    void suffixupdhelper(const szarraytype& ls, T delta, size_t idx = 0) {
        if constexpr (I == dimen) {
            if constexpr (B)
                bit_[idx] = f_(bit_[idx], delta);
            else
                bit_[idx] = g_(bit_[idx], delta);
        } else {
            idx *= szs_[I];
            for (size_t l = ls[I]; l < szs_[I]; l += lsb(l + 1))
                suffixupdhelper<B, I + 1>(ls, delta, idx + l);
        }
    }
    template<bool B, size_t I>
    void rangeupdhelper(
        const szarraytype& ls, const szarraytype& rs, T delta, size_t idx = 0)
    {
        if constexpr (I == dimen) {
            if constexpr (B)
                bit_[idx] = f_(bit_[idx], delta);
            else
                bit_[idx] = g_(bit_[idx], delta);
        } else {
            idx *= szs_[I];
            for (size_t l = ls[I]; l < szs_[I]; l += lsb(l + 1))
                rangeupdhelper<B, I + 1>(ls, rs, delta, idx + l);
            for (size_t r = rs[I]; r < szs_[I]; r += lsb(r + 1))
                rangeupdhelper<!B, I + 1>(ls, rs, delta, idx + r);
        }
    }
    szarraytype szs_;
    vector<T> bit_;
    F f_; G g_;
};

// https://arxiv.org/pdf/1311.6093v4.pdf
template<typename T, T Id, typename F, typename G, typename H, size_t... Idxs>
class ndfwtreerrimpl<T, Id, F, G, H, index_sequence<Idxs...>> {
public:
    static constexpr size_t dimen = sizeof...(Idxs);
    static constexpr size_t k = 1 << dimen;
    static constexpr T identity = Id;
    using valuetype = T;
    using ftype = F;
    using gtype = G;
    using htype = H;
    using szarraytype = array<size_t, dimen>;
    template<typename F1 = F, typename G1 = G, typename H1 = H>
    explicit ndfwtreerrimpl(
        alwaysszt<Idxs>... szs, F1&& f = F(), G1&& g = G(), H1&& h = H())
        : szs_{szs...}, bits_((... * szs))
        , f_(forward<F1>(f)), g_(forward<G1>(g)), h_(forward<H1>(h)) {}
    template<typename F1 = F, typename G1 = G, typename H1 = H>
    explicit ndfwtreerrimpl(
        const szarraytype& szs, F1&& f = F(), G1&& g = G(), H1&& h = H())
        : ndfwtreerrimpl(
            szs, forward<F1>(f), forward<G1>(g), forward<H1>(h),
            make_index_sequence<dimen>()) {}
    template<typename F1 = F, typename G1 = G, typename H1 = H>
    explicit ndfwtreerrimpl(
        const ndvectortype<T, dimen>& v,
        F1&& f = F(), G1&& g = G(), H1&& h = H())
        : ndfwtreerrimpl(
            ndvectorsizes(v), forward<F1>(f), forward<G1>(g), forward<H1>(h))
    {
        szarraytype is;
        preprocess<0>(is, v);
        buildbits<0>(bits_, szs_, f_);
    }
    T prefixqry(alwaysszt<Idxs>... rs) const {
        return prefixqry(array{rs...});
    }
    T prefixqry(const szarraytype& rs) const {
        return prefixqryhelper<0, 0, k>(rs);
    }
    T rangeqry(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs) const {
        return rangeqry(array{ls...}, array{rs...});
    }
    T rangeqry(const szarraytype& ls, const szarraytype& rs) const {
        return rangeqryhelper<0, 0, k>(ls, rs);
    }
    void suffixupd(alwaysszt<Idxs>... ls, T delta) {
        suffixupd(array{ls...}, delta);
    }
    void suffixupd(const szarraytype& ls, T delta) {
        suffixupdhelper<true, 0, 0, k>(ls, delta);
    }
    void suffixupdinv(alwaysszt<Idxs>... ls, T delta) {
        suffixupdinv(array{ls...}, delta);
    }
    void suffixupdinv(const szarraytype& ls, T delta) {
        suffixupdhelper<false, 0, 0, k>(ls, delta);
    }
    void rangeupd(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs, T delta) {
        rangeupd(array{ls...}, array{rs...}, delta);
    }
    void rangeupd(const szarraytype& ls, const szarraytype& rs, T delta) {
        rangeupdhelper<true, 0, 0, k>(ls, rs, delta);
    }
    void rangeupdinv(alwaysszt<Idxs>... ls, alwaysszt<Idxs>... rs, T delta) {
        rangeupdinv(array{ls...}, array{rs...}, delta);
    }
    void rangeupdinv(const szarraytype& ls, const szarraytype& rs, T delta) {
        rangeupdhelper<false, 0, 0, k>(ls, rs, delta);
    }
    template<size_t I = 0>
    size_t size() const { return szs_[I]; }
    F f() const { return f_; }
    G g() const { return g_; }
    H h() const { return h_; }
private:
    template<typename F1, typename G1, typename H1, size_t... Idxs1>
    ndfwtreerrimpl(
        const szarraytype& szs, F1&& f, G1&& g, H1&& h,
        index_sequence<Idxs1...>)
        : ndfwtreerrimpl(
            szs[Idxs1]..., forward<F1>(f), forward<G1>(g), forward<H1>(h)) {}
    template<size_t I, typename U>
    void preprocess(szarraytype& is, const U& u) {
        if constexpr (I == dimen) {
            preprocesshelper<true, 0, 0, k>(is, u);
        } else {
            for (is[I] = 0; is[I] < u.size(); ++is[I])
                preprocess<I + 1>(is, u[is[I]]);
        }
    }
    template<bool B, size_t I, size_t L, size_t R>
    void preprocesshelper(const szarraytype& is, T delta, size_t idx = 0) {
        if constexpr (I == dimen) {
            if constexpr (B)
                bits_[idx][L] = f_(bits_[idx][L], delta);
            else
                bits_[idx][L] = g_(bits_[idx][L], delta);
        } else {
            constexpr size_t mid = (L + R) / 2;
            idx *= szs_[I];
            preprocesshelper<B, I + 1, L, mid>(is, delta, idx + is[I]);
            preprocesshelper<B, I + 1, mid, R>(
                is, h_(delta, is[I]), idx + is[I]);
            if (is[I] + 1 < szs_[I]) {
                preprocesshelper<!B, I + 1, L, mid>(is, delta, idx + is[I] + 1);
                preprocesshelper<!B, I + 1, mid, R>(
                    is, h_(delta, is[I] + 1), idx + is[I] + 1);
            }
        }
    }
    template<size_t I, size_t L, size_t R>
    T prefixqryhelper(const szarraytype& rs, size_t idx = 0) const {
        if constexpr (I == dimen) {
            return bits_[idx][L];
        } else {
            constexpr size_t mid = (L + R) / 2;
            idx *= szs_[I];
            T res = Id;
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                res = f_(res, prefixqryhelper<I + 1, L, mid>(rs, idx + r - 1));
            res = h_(res, rs[I]);
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                res = g_(res, prefixqryhelper<I + 1, mid, R>(rs, idx + r - 1));
            return res;
        }
    }
    template<size_t I, size_t L, size_t R>
    T rangeqryhelper(
        const szarraytype& ls, const szarraytype& rs, size_t idx = 0) const
    {
        if constexpr (I == dimen) {
            return bits_[idx][L];
        } else {
            constexpr size_t mid = (L + R) / 2;
            idx *= szs_[I];
            T rres = Id;
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                rres = f_(
                    rres, rangeqryhelper<I + 1, L, mid>(ls, rs, idx + r - 1));
            rres = h_(rres, rs[I]);
            for (size_t r = rs[I]; r > 0; r -= lsb(r))
                rres = g_(
                    rres, rangeqryhelper<I + 1, mid, R>(ls, rs, idx + r - 1));
            T lres = Id;
            for (size_t l = ls[I]; l > 0; l -= lsb(l))
                lres = f_(
                    lres, rangeqryhelper<I + 1, L, mid>(ls, rs, idx + l - 1));
            lres = h_(lres, ls[I]);
            for (size_t l = ls[I]; l > 0; l -= lsb(l))
                lres = g_(
                    lres, rangeqryhelper<I + 1, mid, R>(ls, rs, idx + l - 1));
            return g_(rres, lres);
        }
    }
    template<bool B, size_t I, size_t L, size_t R>
    void suffixupdhelper(const szarraytype& ls, T delta, size_t idx = 0) {
        if constexpr (I == dimen) {
            if constexpr (B)
                bits_[idx][L] = f_(bits_[idx][L], delta);
            else
                bits_[idx][L] = g_(bits_[idx][L], delta);
        } else {
            constexpr size_t mid = (L + R) / 2;
            idx *= szs_[I];
            for (size_t l = ls[I]; l < szs_[I]; l += lsb(l + 1)) {
                suffixupdhelper<B, I + 1, L, mid>(ls, delta, idx + l);
                suffixupdhelper<B, I + 1, mid, R>(
                    ls, h_(delta, ls[I]), idx + l);
            }
        }
    }
    template<bool B, size_t I, size_t L, size_t R>
    void rangeupdhelper(
        const szarraytype& ls, const szarraytype& rs, T delta, size_t idx = 0)
    {
        if constexpr (I == dimen) {
            if constexpr (B)
                bits_[idx][L] = f_(bits_[idx][L], delta);
            else
                bits_[idx][L] = g_(bits_[idx][L], delta);
        } else {
            constexpr size_t mid = (L + R) / 2;
            idx *= szs_[I];
            for (size_t l = ls[I]; l < szs_[I]; l += lsb(l + 1)) {
                rangeupdhelper<B, I + 1, L, mid>(ls, rs, delta, idx + l);
                rangeupdhelper<B, I + 1, mid, R>(
                    ls, rs, h_(delta, ls[I]), idx + l);
            }
            if (rs[I] < szs_[I]) {
                for (size_t r = rs[I]; r < szs_[I]; r += lsb(r + 1)) {
                    rangeupdhelper<!B, I + 1, L, mid>(ls, rs, delta, idx + r);
                    rangeupdhelper<!B, I + 1, mid, R>(
                        ls, rs, h_(delta, rs[I]), idx + r);
                }
            }
        }
    }
    szarraytype szs_;
    vector<array<T, k>> bits_;
    F f_; G g_; H h_;
};

} // namespace detail

// range-point fenwick tree (range-query + point-update)
template<
    typename T, T Id = 0,
    typename F = plus<T>, typename G = minus<T>,
    size_t N = 1>
class ndfwtree
    : public detail::ndfwtreeimpl<T, Id, F, G, make_index_sequence<N>>
{
    using impl = detail::ndfwtreeimpl<T, Id, F, G, make_index_sequence<N>>;
public:
    using impl::impl;
};

// point-range fenwick tree (point-query + range-update)
template<
    typename T, T Id = 0,
    typename F = plus<T>, typename G = minus<T>,
    size_t N = 1>
class ndfwtreepr
    : public detail::ndfwtreeprimpl<T, Id, F, G, make_index_sequence<N>>
{
    using impl = detail::ndfwtreeprimpl<T, Id, F, G, make_index_sequence<N>>;
public:
    using impl::impl;
};

// range-range fenwick tree (range-query + range-update)
template<
    typename T, T Id = 0,
    typename F = plus<T>, typename G = minus<T>, typename H = multiplies<T>,
    size_t N = 1>
class ndfwtreerr
    : public detail::ndfwtreerrimpl<T, Id, F, G, H, make_index_sequence<N>>
{
    using impl = detail::ndfwtreerrimpl<
        T, Id, F, G, H, make_index_sequence<N>>;
public:
    using impl::impl;
};

// helpful aliases
template<typename T>
struct powfn {
    T operator()(T lhs, size_t rhs) const {
        return static_cast<T>(pow(lhs, rhs));
    }
};
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

int main() {
    {
        vector<uintmax_t> v(100);
        iota(v.begin(), v.end(), 1);
        fwtreemul<uintmax_t> fw(move(v));
        auto choose = [&fw](int n, int r) {
            return fw.rangeqry(n - r, n) / fw.prefixqry(r);
        };
        cout << "20 choose 5 = " << choose(20, 5) << "\n";
    }
    cout << "\n";
    {
        auto xorntimes = [](int x, size_t n) { return n & 1 ? x : 0; };
        ndfwtreerr<int, 0, bit_xor<int>, bit_xor<int>, decltype(xorntimes)>
            fw(10, bit_xor<int>(), bit_xor<int>(), xorntimes);
        fw.rangeupd(4, 7, 7);
        fw.rangeupd(2, 6, 7);
        for (size_t i = 0; i < fw.size(); ++i)
            cout << fw.rangeqry(i, i + 1) << " ";
        cout << "\n";
    }
    cout << "\n";
    {
        fwtreepr2d<int> fw(10, 10);
        fw.rangeupd(0, 0, 10, 10, 5);
        fw.rangeupdinv(2, 1, 5, 4, 5);
        fw.rangeupdinv(2, 6, 5, 9, 5);
        fw.rangeupdinv(7, 2, 9, 8, 5);
        for (size_t i = 0; i < fw.size(); ++i) {
            for (size_t j = 0; j < fw.size<1>(); ++j)
                cout << fw.pointqry(i, j) << " ";
            cout << "\n";
        }
    }
    cout << "\n";
    {
        fwtreerr2d<int> fw({
            {1,  1,  0,  1, 1},
            {1,  0, -1,  0, 1},
            {0, -1, -2, -1, 0},
            {1,  0, -1,  0, 1},
            {1,  1,  0,  1, 1}
        });
        fw.rangeupd(1, 1, 4, 4, 1);
        fw.rangeupd(2, 0, 3, 5, 1);
        fw.rangeupd(0, 2, 5, 3, 1);
        for (size_t i = 0; i < fw.size(); ++i) {
            for (size_t j = 0; j < fw.size<1>(); ++j)
                cout << fw.rangeqry(i, j, i + 1, j + 1) << " ";
            cout << "\n";
        }
    }
}
