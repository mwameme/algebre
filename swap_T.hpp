#pragma once

template<class T, class enable = void> class swap_T {
public:
	static inline void swap(T& x, T& y) {
		std::swap(x, y);
		return;
	};
};


template<class T> class swap_T<T, std::void_t<decltype(swap(T(), T()))> > {
public:
	static inline void swap(T& x, T& y) {
		swap(x, y);
		return;
	};
};

template<class T> inline void swap_F(T& x, T& y) {
	swap_T<T>::swap(x, y);
	return;
};
