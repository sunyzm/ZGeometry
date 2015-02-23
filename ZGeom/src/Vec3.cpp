#include "Vec3.h"

namespace ZGeom {

Matrix3::Matrix3()
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            data[i][j] = 0.0;
}

Matrix3::Matrix3(const Matrix3& mat)
{
    std::memcpy((void*)data, (void*)mat.data, sizeof(data));
}

Matrix3& Matrix3::operator=(const Matrix3& mat)
{
    std::memcpy((void*)data, (void*)mat.data, sizeof(data));
    return *this;
}

double& Matrix3::operator()(int i, int j)
{
    assert(0 <= i && i <= 2 && 0 <= j && j <= 2);
    return data[i][j];
}

Matrix3& Matrix3::operator*=(double c)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            data[i][j] *= c;
    return *this;
}

Matrix3& Matrix3::operator+=(const Matrix3& mat1)
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            data[i][j] += mat1.data[i][j];
        }
    }
    return *this;
}

Matrix3 operator+(const Matrix3& mat1, const Matrix3& mat2)
{
    Matrix3 mat;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            mat.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
    return mat;
}

Matrix3 operator-(const Matrix3& mat1, const Matrix3& mat2)
{
    Matrix3 mat;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            mat.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
    return mat;
}

Vec3d operator*(const Vec3d& vec, const Matrix3& mat1)
{
    Vec3d ret;
    Matrix3 mat(mat1);

    ret.x = vec.x * mat(0, 0) + vec.y * mat(1, 0) + vec.z * mat(2, 0);
    ret.y = vec.x * mat(0, 1) + vec.y * mat(1, 1) + vec.z * mat(2, 1);
    ret.z = vec.x * mat(0, 2) + vec.y * mat(1, 2) + vec.z * mat(2, 2);

    return ret;
}

Vec3d operator*(const Matrix3& mat1, const Vec3d& vec)
{
    Vec3d ret;
    Matrix3 mat(mat1);

    ret.x = mat(0, 0) * vec.x + mat(0, 1) * vec.y + mat(0, 2) * vec.z;
    ret.x = mat(1, 0) * vec.x + mat(1, 1) * vec.y + mat(1, 2) * vec.z;
    ret.x = mat(2, 0) * vec.x + mat(2, 1) * vec.y + mat(2, 2) * vec.z;
    return ret;
}

Matrix3 vector3DMultiply(const Vec3d& v1, const Vec3d& v2)
{
    Matrix3 mat;

    mat.data[0][0] = v1.x * v2.x;
    mat.data[0][1] = v1.x * v2.y;
    mat.data[0][2] = v1.x * v2.z;
    mat.data[1][0] = v1.y * v2.x;
    mat.data[1][1] = v1.y * v2.y;
    mat.data[1][2] = v1.y * v2.z;
    mat.data[2][0] = v1.z * v2.x;
    mat.data[2][1] = v1.z * v2.y;
    mat.data[2][2] = v1.z * v2.z;

    return mat;
}


}   // end of namespace

