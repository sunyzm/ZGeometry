#pragma once
#include <vector>
#include <map>
#include <string>

#define ID_PROPERTY_NOT_SPECIFIED -1

class MeshProcessor;

class MeshProperty
{
public:
	friend class MeshProcessor;

	MeshProperty() : id(ID_PROPERTY_NOT_SPECIFIED), name("unnamed") {}
	void setIDandName(int newID, const std::string& newName ) { id = newID; name = newName; }
protected:
	int id;
	std::string name;
};

class MeshHandle
{
public:
	int index;
	Vector3D position;
	MeshHandle(int i, double vx, double vy, double vz) : index(i) { position = Vector3D(vx, vy, vz); }
	MeshHandle(int i, const Vector3D& pos) : index(i), position(pos) {}
};

class MeshFeature
{
public:
	int index;
	int scale;
	double value;
	MeshFeature(int i, int s) : index(i), scale(s) {}
};

class MeshFeatureList : public MeshProperty
{
public:
	int featureType;
	std::vector<MeshFeature> m_vFeatures;
};

class MeshFunction : public MeshProperty
{
	int m_size;		// vertex number of the mesh
public:
	MeshFunction() { m_size = 0; }
	MeshFunction(int s) { m_size = 0; m_function.resize(s); }
	void setSize(int s) { m_size = s; m_function.resize(s); }
	const std::vector<double>& getMeshFunction() { return m_function; }
	static double InnerProduct(const MeshFunction& f1, const MeshFunction& f2);
	double norm() const;

	double& operator[](int idx);
protected:
	std::vector<double> m_function;
};

class SparseManifoldFunction : public MeshProperty
{
public:
	std::vector<int> m_funcIdx;
	std::vector<double> m_funcVal;
};

class MeshKernel : public MeshProperty
{

};

class MeshProcessor
{
public:
	MeshProcessor();
	void setMesh(CMesh* newMesh);
	bool addProperty(MeshProperty* newProperty);
	MeshProperty* retrievePropertyByID(int rid);
	MeshProperty* retrievePropertyByName(const std::string& rn);
	void removePropertyByID(int rid);
	void removePropertyByName(const std::string& rn);
protected:
	CMesh* mesh;
	int m_size;

	std::vector<MeshProperty*> vProperties;
};