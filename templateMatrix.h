#ifndef TEMPL_MATRICES_HH
#define TEMPL_MATRICES_HH


template<typename _T>
inline void ez_swap(_T& A, _T& B) {
    _T tmp = A;
    A = B;
    B = tmp;
}


template<typename _T, int ROWS, int COLS>//, uint ROWS, uint COLS>
class Matrix {
private:
public:
#ifdef NAME_FOLLOW
    std::string name;
#endif
    _T buffer[ROWS * COLS];
    int row_stride;
    int col_stride;

    void fast_row_update(uint row_number, const _T* src) {
        memcpy(buffer + row_number * row_stride, src, COLS * sizeof(_T));
    }
    inline const _T* get_const_elem(const uint r, const int c)const {
        return &buffer[r * row_stride + c * col_stride];
    }
    inline _T* get_elem(const uint r, const int c) {
        return &buffer[r * row_stride + c * col_stride];
    }
    void row_update(const uint row_number, const _T* const src) {
        for (int i = 0;i < COLS;i++) {
            //buffer[row_number * row_stride + i * col_stride] = src[i];
            *get_elem(row_number, i) = src[i];
        }
    }
    void col_update(const uint col_number, const _T* const src) {
        for (int i = 0;i < ROWS;i++) {
            //buffer[row_number * row_stride + i * col_stride] = src[i];
            *get_elem(i, col_number) = src[i];
        }
    }


    Matrix() {
        col_stride = 1;
        row_stride = COLS;
        memset(buffer, 0, ROWS * COLS * sizeof(_T));
    }
#ifdef NAME_FOLLOW
    void setname(const std::string& name) {
        this->name = name;
    }
#endif


    Matrix(const std::vector<std::vector<_T>>& elems) {


        col_stride = 1;
        row_stride = COLS;

        for (int i = 0;i < ROWS;i++) {
            fast_row_update(i, elems[i].data());
        }
    }

    void inplace_col_mul(const int col, const _T k) {
        for (int r = 0;r < ROWS;r++) {
            (*get_elem(r, col)) *= k;
        }
    }
    void inplace_row_mul(const int row, const _T k) {
        for (int c = 0;c < COLS;c++) {
            (*get_elem(row, c)) *= k;
        }
    }
    void weighted_rows(int dst_row, int row1, _T w1, int row2, _T w2) {
        for (int c = 0;c < COLS;c++) { //columns
            *get_elem(dst_row, c) = w1 * *get_elem(row1, c) + w2 * *get_elem(row2, c);
        }
    }
    void weighted_cols(int dst_col, int col1, _T w1, int col2, _T w2) {
        for (int r = 0;r < ROWS;r++) {
            *get_elem(r, dst_col) = w1 * *get_elem(r, col1) + w2 * *get_elem(r, col2);
        }
    }
    void add_weighted_row(const int dst_row, const int row, const _T w) {
        for (int c = 0;c < COLS;c++) {
            *get_elem(dst_row, c) += w * *get_elem(row, c);
        }
    }
    void add_weighted_col(const int dst_col, const int col, const _T w) {
        for (int r = 0;r < ROWS;r++) {
            *get_elem(r, dst_col) += w * *get_elem(r, col);
        }
    }

    ~Matrix() {
#ifdef NAME_FOLLOW
        std::cout << "Freed " << name << "\n";
#endif
    }


};
template<int ROWS, int COLS>
void Identity(Matrix<float, ROWS, COLS>& M) {
    memset(M.buffer, 0, ROWS * COLS * sizeof(float));
    for (int i = 0;i < std::min(ROWS, COLS);i++) {
        *M.get_elem(i, i) = 1.f;
    }
}

template<int X, int Y, int Z>
Matrix<float, X, Z> operator*(const Matrix<float, X, Y>& A, const Matrix<float, Y, Z>& B) {

    Matrix<float, X, Z> tmp;
#ifdef NAME_FOLLOW
    tmp.setname("(" + A.name + "*" + B.name + ")");
#endif
    for (int r = 0;r < X;r++) {
        for (int c = 0;c < Z;c++) {
            //Row * Col
            float& current_element = *tmp.get_elem(r, c);
            //current_element = 0;
            for (int k = 0;k < Y;k++) {
                current_element += *A.get_const_elem(r, k) * (*B.get_const_elem(k, c));
            }

        }
    }

    return tmp;
}

template<int X, int Y>
Matrix<float, X, Y> operator+(const Matrix<float, X, Y>& A, const Matrix<float, X, Y>& B) {

    Matrix<float, X, Y> tmp;
#ifdef NAME_FOLLOW
    tmp.setname("(" + A.name + "+" + B.name + ")");
#endif
    for (int r = 0;r < X;r++) {
        for (int c = 0;c < Y;c++) {
            //Row * Col
            float& current_element = *tmp.get_elem(r, c);
            //current_element = 0;
            current_element = *A.get_const_elem(r, c) + (*B.get_const_elem(r, c));
        }
    }

    return tmp;
}

template<int X, int Y>
Matrix<float, X, Y> operator-(const Matrix<float, X, Y>& A, const Matrix<float, X, Y>& B) {

    Matrix<float, X, Y> tmp;
#ifdef NAME_FOLLOW
    tmp.setname("(" + A.name + "-" + B.name + ")");
#endif
    for (int r = 0;r < X;r++) {
        for (int c = 0;c < Y;c++) {
            //Row * Col
            float& current_element = *tmp.get_elem(r, c);
            //current_element = 0;
            current_element = *A.get_const_elem(r, c) - (*B.get_const_elem(r, c));
        }
    }

    return tmp;
}

template<int COLS, int ROWS>
Matrix<float, COLS, ROWS> Transpose(const Matrix<float, ROWS, COLS>& M) {
    Matrix<float, COLS, ROWS> tmp(M);
    tmp.row_stride = 1;
    tmp.col_stride = ROWS;
    return tmp;
}

template<int N>
Matrix<float, N, N> Inverse(const Matrix<float, N, N>& M) {
    // tmp.setname("tmp");
    Matrix<float, N, N> Inv;

    for (int i = 0;i < N;i++) {
        *Inv.get_elem(i, i) = 1;
    }

    Matrix<float, N, N> tmp = M;

#ifdef NAME_FOLLOW
    tmp.setname("TMP");
#endif
    for (int c = 0;c < N - 1;c++) {
        for (int r = c + 1;r < N;r++) {
            float kc = -(*tmp.get_elem(r, c)) / (*tmp.get_elem(c, c));
            tmp.add_weighted_row(r, c, kc);
            Inv.add_weighted_row(r, c, kc);
        }
    }

    for (int c = N - 1;c > 0;c--) {
        for (int r = c - 1;r >= 0;r--) {
            float kc = -(*tmp.get_elem(r, c)) / (*tmp.get_elem(c, c));
            tmp.add_weighted_row(r, c, kc);
            Inv.add_weighted_row(r, c, kc);

        }
    }

    for (int c = 0;c < N;c++) {
        Inv.inplace_row_mul(c, float(1) / (*tmp.get_elem(c, c)));
    }
    //Inv.Print();
    return Inv;

}

#endif
