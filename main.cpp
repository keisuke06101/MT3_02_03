#include <Novice.h>
#include <cmath>
#include <math.h>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <cassert>
#include <numbers>
#include <ImGuiManager.h>

const char kWindowTitle[] = "GC2A_16_タナカケイスケ";

struct Matrix4x4
{
	float m[4][4];
};


struct Vector3
{
	float x;
	float y;
	float z;
};

struct Line {
	Vector3 origin;  //!< 始点
	Vector3 dift;    //!< 終点への差分ベクトル

	static constexpr float KTMin = std::numeric_limits<float>::lowest();
	static constexpr float KTMax = (std::numeric_limits<float>::max)();
};

struct Ray {
	Vector3 origin;  //!< 始点
	Vector3 dift;    //!< 終点への差分ベクトル
};

struct Segment {
	Vector3 origin;  //!< 始点
	Vector3 diff;    //!< 終点への差分ベクトル
	uint32_t color;

	static constexpr float KTMin = 0.0f;
	static constexpr float KTMax = 1.0f;
};

struct Sphere
{
	Vector3 center;
	float radius;
	uint32_t color;
};

struct Plane {
	Vector3 normal;
	float distance;
};

struct Triangle {
	Vector3 vertices[3]; //!< 頂点
};

// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m)
{
	Matrix4x4 resultInverse;
	float determinant = m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2] -

		m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2] -

		m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2] +

		m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2] +

		m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2] -

		m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2] -

		m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0] +

		m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];


	assert(determinant != 0.0f);
	float determinantRecp = 1.0f / determinant;

	//1行目
	resultInverse.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[1][3] * m.m[2][1] * m.m[3][2] -
		m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[1][1] * m.m[2][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[0][1] = -(m.m[0][1] * m.m[2][2] * m.m[3][3] +
		m.m[0][2] * m.m[2][3] * m.m[3][1] +
		m.m[0][3] * m.m[2][1] * m.m[3][2] -
		m.m[0][3] * m.m[2][2] * m.m[3][1] -
		m.m[0][2] * m.m[2][1] * m.m[3][3] -
		m.m[0][1] * m.m[2][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[3][2] -
		m.m[0][3] * m.m[1][2] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[0][3] = -(m.m[0][1] * m.m[1][2] * m.m[2][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][1] +
		m.m[0][3] * m.m[1][1] * m.m[2][2] -
		m.m[0][3] * m.m[1][2] * m.m[2][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][3] -
		m.m[0][1] * m.m[1][3] * m.m[2][2]) * determinantRecp;

	//2行目
	resultInverse.m[1][0] = -(m.m[1][0] * m.m[2][2] * m.m[3][3] +
		m.m[1][2] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][2] -
		m.m[1][3] * m.m[2][2] * m.m[3][0] -
		m.m[1][2] * m.m[2][0] * m.m[3][3] -
		m.m[1][0] * m.m[2][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] +
		m.m[0][2] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][2] -
		m.m[0][3] * m.m[2][2] * m.m[3][0] -
		m.m[0][2] * m.m[2][0] * m.m[3][3] -
		m.m[0][0] * m.m[2][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[1][2] = -(m.m[0][0] * m.m[1][2] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][2] -
		m.m[0][3] * m.m[1][2] * m.m[3][0] -
		m.m[0][2] * m.m[1][0] * m.m[3][3] -
		m.m[0][0] * m.m[1][3] * m.m[3][2]) * determinantRecp;

	resultInverse.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] -
		m.m[0][3] * m.m[1][2] * m.m[2][0] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] -
		m.m[0][0] * m.m[1][3] * m.m[2][2]) * determinantRecp;

	//3行目
	resultInverse.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] +
		m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][1] -
		m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[1][0] * m.m[2][3] * m.m[3][1]) * determinantRecp;

	resultInverse.m[2][1] = -(m.m[0][0] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][1] -
		m.m[0][3] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][3] -
		m.m[0][0] * m.m[2][3] * m.m[3][1]) * determinantRecp;

	resultInverse.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] +
		m.m[0][1] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][1] -
		m.m[0][3] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][3] -
		m.m[0][0] * m.m[1][3] * m.m[3][1]) * determinantRecp;

	resultInverse.m[2][3] = -(m.m[0][0] * m.m[1][1] * m.m[2][3] +
		m.m[0][1] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][1] -
		m.m[0][3] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][3] -
		m.m[0][0] * m.m[1][3] * m.m[2][1]) * determinantRecp;

	//4行目
	resultInverse.m[3][0] = -(m.m[1][0] * m.m[2][1] * m.m[3][2] +
		m.m[1][1] * m.m[2][2] * m.m[3][0] +
		m.m[1][2] * m.m[2][0] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][2] -
		m.m[1][0] * m.m[2][2] * m.m[3][1]) * determinantRecp;

	resultInverse.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] +
		m.m[0][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][2] * m.m[2][0] * m.m[3][1] -
		m.m[0][2] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][2] -
		m.m[0][0] * m.m[2][2] * m.m[3][1]) * determinantRecp;

	resultInverse.m[3][2] = -(m.m[0][0] * m.m[1][1] * m.m[3][2] +
		m.m[0][1] * m.m[1][2] * m.m[3][0] +
		m.m[0][2] * m.m[1][0] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][2] -
		m.m[0][0] * m.m[1][2] * m.m[3][1]) * determinantRecp;

	resultInverse.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] +
		m.m[0][1] * m.m[1][2] * m.m[2][0] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][2] -
		m.m[0][0] * m.m[1][2] * m.m[2][1]) * determinantRecp;

	return resultInverse;
}

// 透視射影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
	Matrix4x4 result;
	float cot = 1.0f / std::tan(fovY / 2.0f);

	result = { (1.0f / aspectRatio) * cot, 0.0f,    0.0f,                                       0.0f,
				0.0f,                      cot,     0.0f,                                       0.0f,
				0.0f,                      0.0f,    farClip / (farClip - nearClip),             1.0f,
				0.0f,                      0.0f,    -nearClip * farClip / (farClip - nearClip), 0.0f };

	return result;
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
	Matrix4x4 result;
	result = { width / 2,        0.0f,           0.0f,              0.0f,
			   0.0f,             -height / 2,      0.0f,                0.0f,
			   0.0f,             0.0f,             maxDepth - minDepth, 0.0f,
			   left + width / 2, top + height / 2, minDepth,            1.0f };

	return result;
}


// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix)
{
	Vector3 result;
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;
	return result;
}


//スカラー倍
Vector3 Multiply(float scalar, const Vector3& v)
{
	return{ scalar * v.x, scalar * v.y, scalar * v.z };
};

//減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2)
{
	return{ v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
};

//加算
Vector3 Add(const Vector3& v1, const Vector3& v2)
{
	return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
};
//長さ（ノルム）
float Length(const Vector3& v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
};
//正規化
Vector3 Normalize(const Vector3& v)
{
	float length = Length(v);
	return { v.x / length, v.y / length, v.z / length };
};

//内積
float Dot(const Vector3& v1, const Vector3& v2)
{
	return { v1.x * v2.x + v1.y * v2.y + v1.z * v2.z };
};

// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
	Matrix4x4 resultMultiply;
	resultMultiply.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
	resultMultiply.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
	resultMultiply.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
	resultMultiply.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

	resultMultiply.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
	resultMultiply.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
	resultMultiply.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
	resultMultiply.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

	resultMultiply.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
	resultMultiply.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
	resultMultiply.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
	resultMultiply.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

	resultMultiply.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
	resultMultiply.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
	resultMultiply.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
	resultMultiply.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];
	return resultMultiply;
}



/// <summary>
/// 拡大縮小行列
/// </summary>
/// <param name="scale"></param>
/// <returns></returns>
Matrix4x4 Scale(const Vector3& scale) {
	Matrix4x4 result;
	result = { scale.x, 0.0f,    0.0f,    0.0f,
			  0.0f,    scale.y, 0.0f,    0.0f,
			  0.0f,    0.0f,    scale.z, 0.0f,
			  0.0f,    0.0f,    0.0f,    1.0f };
	return result;
}


// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian) {
	float sin = std::sin(radian);
	float cos = std::cos(radian);

	Matrix4x4 result{ 1.0f, 0.0f, 0.0f, 0.0f,
					  0.0f, cos,  sin,  0.0f,
					  0.0f, -sin, cos,  0.0f,
					  0.0f, 0.0f, 0.0f, 1.0f };
	return result;
}

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	float sin = std::sin(radian);
	float cos = std::cos(radian);

	Matrix4x4 result{ cos,  0.0f, -sin, 0.0f,
					  0.0f, 1.0f, 0.0f, 0.0f,
					  sin,  0.0f, cos,  0.0f,
					  0.0f, 0.0f, 0.0f, 1.0f };
	return result;
}

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	float sin = std::sin(radian);
	float cos = std::cos(radian);

	Matrix4x4 result{ cos,  sin,  0.0f, 0.0f,
					  -sin, cos,  0.0f, 0.0f,
					  0.0f, 0.0f, 1.0f, 0.0f,
					  0.0f, 0.0f, 0.0f, 1.0f };
	return result;
}

//平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& trans) {
	Matrix4x4 result{ 1.0f,    0.0f,    0.0f,    0.0f,
					  0.0f,    1.0f,    0.0f,    0.0f,
					  0.0f,    0.0f,    1.0f,    0.0f,
					  trans.x, trans.y, trans.z, 1.0f };
	return result;
}


// 代入演算子
Matrix4x4& operator*=(Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			for (size_t k = 0; k < 4; k++) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	m1 = result;
	return m1;
}

// 二項演算子
Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2)
{
	Matrix4x4 result = m1;

	return result *= m2;
}

// 3次元アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	// スケーリング
	Matrix4x4 matScale = Scale(scale);

	Matrix4x4 matRotateX = MakeRotateXMatrix(rotate.x);
	Matrix4x4 matRotateY = MakeRotateYMatrix(rotate.y);
	Matrix4x4 matRotateZ = MakeRotateZMatrix(rotate.z);

	Matrix4x4 matRot = matRotateX * matRotateY * matRotateZ;

	Matrix4x4 matTranslate = MakeTranslateMatrix(translate);

	Matrix4x4 result = matScale * matRot * matTranslate;

	return result;
}

// クロス変換
Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
	Vector3 result;

	result = { (v1.y * v2.z - v1.z * v2.y) ,(v1.z * v2.x - v1.x * v2.z) , (v1.x * v2.y - v1.y * v2.x) };

	return result;

}

// グリッド表示
void DrawGrid(const Matrix4x4& viewProjectMatrix, const Matrix4x4& viewportMatrix)
{
	const float kGridHalfWidth = 2.0f;    // Gridの半分の幅

	const uint32_t kSubdivision = 10;     // 分割数

	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);  // 1つ分の長さ

	// 奥から手前への線を引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex)
	{
		float x = -kGridHalfWidth + (xIndex * kGridEvery);

		Vector3 stert = { x,0.0f,-kGridHalfWidth };

		Vector3 stertScr = Transform(Transform(stert, viewProjectMatrix), viewportMatrix);

		Vector3 end = { x,0.0f, kGridHalfWidth };

		Vector3 endScr = Transform(Transform(end, viewProjectMatrix), viewportMatrix);

		Novice::DrawLine(int(stertScr.x), int(stertScr.y), int(endScr.x), int(endScr.y), x == 0.0f ? BLACK : 0xAAAAAAFF);

	}

	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex)
	{
		float z = -kGridHalfWidth + (zIndex * kGridEvery);

		Vector3 stert = { -kGridHalfWidth,0.0f,z };

		Vector3 stertScr = Transform(Transform(stert, viewProjectMatrix), viewportMatrix);

		Vector3 end = { kGridHalfWidth,0.0f,z };

		Vector3 endScr = Transform(Transform(end, viewProjectMatrix), viewportMatrix);

		Novice::DrawLine(int(stertScr.x), int(stertScr.y), int(endScr.x), int(endScr.y), z == 0.0f ? BLACK : 0xAAAAAAFF);
	}

}

//スフィア表示
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	float pi = std::numbers::pi_v<float>;

	const uint32_t kSubdivision = 12;  // 分割数

	const float kLonEvery = pi * 2.0f / float(kSubdivision);        // 緯度分割1つ分の角度

	const float kLatEvery = pi / float(kSubdivision);               // 経度分割1つ分の角度

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex)
	{
		float lat = -pi / 2.0f + kLatEvery * latIndex;              // 現在の緯度

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex)
		{
			float lon = lonIndex * kLonEvery;                         // 現在の経度

			Vector3 a = {
				sphere.center.x + sphere.radius * std::cos(lat) * std::cos(lon),
				sphere.center.y + sphere.radius * std::sin(lat),
				sphere.center.z + sphere.radius * std::cos(lat) * std::sin(lon) };

			Vector3 b = {
				sphere.center.x + sphere.radius * std::cos(lat + kLatEvery) * std::cos(lon),
				sphere.center.y + sphere.radius * std::sin(lat + kLatEvery),
				sphere.center.z + sphere.radius * std::cos(lat + kLatEvery) * std::sin(lon) };

			Vector3 c = {
				sphere.center.x + sphere.radius * std::cos(lat) * std::cos(lon + kLonEvery),
				sphere.center.y + sphere.radius * std::sin(lat),
				sphere.center.z + sphere.radius * std::cos(lat) * std::sin(lon + kLonEvery) };


			// 線を描く
			Vector3 scrA = Transform(Transform(a, viewProjectMatrix), viewportMatrix);
			Vector3 scrB = Transform(Transform(b, viewProjectMatrix), viewportMatrix);
			Vector3 scrC = Transform(Transform(c, viewProjectMatrix), viewportMatrix);

			Novice::DrawLine(int(scrA.x), int(scrA.y), int(scrB.x), int(scrB.y), color);
			Novice::DrawLine(int(scrB.x), int(scrB.y), int(scrC.x), int(scrC.y), color);
		}
	}
}

// 平面
Vector3 Perpendicular(const Vector3& vector)
{
	if (vector.x != 0.0f || vector.y != 0.0f)
	{
		return { -vector.y, vector.x, 0.0f };
	}
	return { 0.0f, -vector.z, vector.y };
}

// 平面の描画
void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4]{};
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	perpendiculars[1] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	perpendiculars[3] = { -perpendiculars[2].x, -perpendiculars[2].y, -perpendiculars[2].z };
	Vector3 points[4]{};

	for (int32_t index = 0; index < 4; ++index)
	{
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}

	Novice::DrawLine(int(points[0].x), int(points[0].y), int(points[2].x), int(points[2].y), color);
	Novice::DrawLine(int(points[1].x), int(points[1].y), int(points[3].x), int(points[3].y), color);
	Novice::DrawLine(int(points[2].x), int(points[2].y), int(points[1].x), int(points[1].y), color);
	Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[0].x), int(points[0].y), color);
}

Vector3 Project(const Vector3& v1, const Vector3& v2)
{
	float v25qLegth = Dot(v2, v2);
	float dot = Dot(v1, v2);
	return Multiply(dot / v25qLegth, v2);
}

Vector3 ClosestPoint(const Vector3& point, const Segment& segment)
{
	Vector3 v = Subtract(point, segment.origin);
	float t = Dot(v, segment.diff) / Dot(segment.diff, segment.diff);
	t = std::clamp(t, 0.0f, 1.0f);
	return Add(segment.origin, Multiply(t, segment.diff));
}

void DrawSegment(const Segment& segment, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
	Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);
	Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), color);
}

// 衝突判定処理
bool IsCollision(const Segment& s1, const Plane s2)
{
	// まず垂直判定を行うために、法線と線の内積を求める
	float dot = Dot(s2.normal, s1.diff);

	if (dot != 0.0f)
	{
		// tを求める
		float t = (s2.distance - Dot(s1.origin, s2.normal)) / dot;
		return (Segment::KTMin <= t) && (t <= Segment::KTMax);
	}
	
	return false;
}

bool IsCollisionT(const Triangle& s1, const Segment& s2)
{
	Vector3 v01 = Subtract(s1.vertices[1], s1.vertices[0]);
	Vector3 v12 = Subtract(s1.vertices[2], s1.vertices[1]);
	Vector3 normal = Normalize(Cross(v01, v12));
	Plane plane{ .normal = normal, .distance = Dot(s1.vertices[0], normal) };
	float dot = Dot(plane.normal, s2.diff);
	if (dot == 0.0f)
	{
		return false;
	}
	// tを求める
	float t = (plane.distance - Dot(s2.origin, plane.normal)) / dot;
	if ((Segment::KTMin > t) || (t > Segment::KTMax)) {
		return false;
	}
	Vector3 intersect = Add(s2.origin, Multiply(t, s2.diff));
	Vector3 v1p = Subtract(intersect, s1.vertices[1]);
	if (Dot(Cross(v01, v1p), normal) < 0.0f) {
		return false;
	}

	Vector3 v2p = Subtract(intersect, s1.vertices[2]);
	if (Dot(Cross(v12, v2p), normal) < 0.0f) {
		return false;
	}

	Vector3 v0p = Subtract(intersect, s1.vertices[0]);
	Vector3 v20 = Subtract(s1.vertices[0], s1.vertices[2]);
	if (Dot(Cross(v20, v0p), normal) < 0.0f) {
		return false;
	}

	return true;
}

void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Vector3 screenVertices[3] = {
		Transform(Transform(triangle.vertices[0], viewProjectionMatrix), viewportMatrix),
		Transform(Transform(triangle.vertices[1], viewProjectionMatrix), viewportMatrix),
		Transform(Transform(triangle.vertices[2], viewProjectionMatrix), viewportMatrix),
	};
	Novice::DrawTriangle(
		int(screenVertices[0].x), int(screenVertices[0].y), int(screenVertices[1].x),
		int(screenVertices[1].y), int(screenVertices[2].x), int(screenVertices[2].y), color, kFillModeWireFrame);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Segment segment{ {0.0f, 0.0f, 0.0f}, {-1.0f, 1.0f, 0.0f}, WHITE };
	Vector3 point{ -1.5f, 0.6f, 0.6f };

	Vector3 project = Project(Subtract(point, segment.origin), segment.diff);
	Vector3 closestPoint = ClosestPoint(point, segment);

	Sphere pointSphere{ point, 0.01f };
	Sphere closestPointSphere{ closestPoint, 0.01f };

	Vector3 rotate{};

	Vector3 translate{};

	Vector3 cameraTranslate = { 0.0f,1.9f,-6.49f };

	Vector3 cameraRotate = { 0.26f,0.0f,0.0f };

	Sphere sphere =
	{
		{0.0f,0.0f,0.0f},
		0.5f,
		WHITE
	};

	Plane plane =
	{
		{0.0f,1.0f,0.0f},
		1.0f,
	};

	Triangle triangle =
	{
		{
			{-0.5f, -0.5f, 0.0f},
			{ 0.0f,  0.5f, 0.0f},
			{ 0.5f, -0.5f, 0.0f}
		}
	};

	int kWindowWidth = 1280;

	int kWindowHeight = 720;

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		
		//ImGui
		ImGui::SetNextWindowPos({ 0,0 });
		ImGui::SetNextWindowSize({ 300, 200 });

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("SegmentOrigin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("SegmentDift", &segment.diff.x, 0.01f);
		ImGui::DragFloat3("TriangleV0", &triangle.vertices[0].x, 0.01f);
		ImGui::DragFloat3("TriangleV1", &triangle.vertices[1].x, 0.01f);
		ImGui::DragFloat3("TriangleV2", &triangle.vertices[2].x, 0.01f);
		plane.normal = Normalize(plane.normal);

		ImGui::End();

		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f, 0.0f,0.0f }, translate);

		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);

		Matrix4x4 viewMatrix = Inverse(cameraMatrix);

		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);

		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));

		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		
		Vector3 start = Transform(Transform(segment.origin, worldViewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), worldViewProjectionMatrix), viewportMatrix);

		segment.color = IsCollisionT(triangle, segment) ? RED : WHITE;

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///


		DrawGrid(worldViewProjectionMatrix, viewportMatrix);
		DrawSegment(segment, worldViewProjectionMatrix, viewportMatrix, segment.color);
		DrawTriangle(triangle, worldViewProjectionMatrix, viewportMatrix, WHITE);
		//DrawPlane(plane, worldViewProjectionMatrix, viewportMatrix, WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
