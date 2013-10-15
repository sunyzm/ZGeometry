#pragma once

#include <vector>
#include <map>
#include <string>
#include "Geometry.h"
#include "Mesh.h"

#define ID_PROPERTY_NOT_SPECIFIED 0

class MeshProcessor;

class MeshProperty
{
public:
	friend class MeshProcessor;
	MeshProperty() : id(ID_PROPERTY_NOT_SPECIFIED), name("unnamed") {}
	virtual ~MeshProperty() {}
	void setIDandName(int newID, const std::string& newName ) { id = newID; name = newName; }

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
	int m_index;
	int m_scale;
	int m_note;
	MeshFeature() : m_index(-1), m_scale(-1), m_note(0) {}
	MeshFeature(int i) : m_index(i), m_scale(0), m_note(0) {}
	MeshFeature(int i, int s) : m_index(i), m_scale(s), m_note(0) {}
	virtual ~MeshFeature() {}
};

class MeshFeatureList : public MeshProperty
{
public:
	void addFeature(MeshFeature* mf) { m_vFeatures.push_back(mf); }
	std::vector<MeshFeature*>* getFeatureVector() { return &m_vFeatures; }
	~MeshFeatureList();

	int featureType;
	std::vector<MeshFeature*> m_vFeatures;
};

class MeshFunction : public MeshProperty
{
public:
	MeshFunction() { m_size = 0; }
	MeshFunction(int s) { m_size = s; m_function.resize(s); }
	MeshFunction(int id, const std::string& newName, int size) { setIDandName(id, newName); m_size = size; m_function.resize(size); }
	void setSize(int s) { m_size = s; m_function.resize(s); }
	const std::vector<double>& getMeshFunction_const() { return m_function; }
	std::vector<double>& getMeshFunction() { return m_function; }
	static double InnerProduct(const MeshFunction& f1, const MeshFunction& f2);
	double norm() const;
	void setValue(int idx, double val) { m_function[idx] = val; }
	double& operator[] (int idx);
	double  operator[] (int idx) const;
	void copyValues(const std::vector<double>& values);

protected:
	int m_size;		// vertex number of the mesh
	std::vector<double> m_function;
};

class SparseMeshFunction : public MeshProperty
{
public:
	std::vector<int> m_funcIdx;
	std::vector<double> m_funcVal;
};

class MeshKernel : public MeshProperty
{
protected:
	int m_size;
	int m_nonzeros;
	std::vector<int> m_II;
	std::vector<int> m_JJ;
	std::vector<int> m_VV;
};

class MeshProcessor
{
public:
	MeshProcessor();
	~MeshProcessor();
	void setMesh(CMesh* newMesh);
	std::vector<MeshProperty*> properties() { return vProperties; }
	bool addProperty(MeshProperty* newProperty);
	void replaceProperty(MeshProperty* newProperty);
	MeshProperty* retrievePropertyByID(int rid);
	const MeshProperty* retrievePropertyByID(int rid) const;
	MeshProperty* retrievePropertyByName(const std::string& rn);
	void removePropertyByID(int rid);
	void removePropertyByName(const std::string& rn);
	const CMesh* getMesh_const() const { return mesh; }
	const CMesh* getOriMesh_const() const { return ori_mesh; }
	CMesh* getMesh() const { return mesh; }

protected:
	CMesh* mesh;
	CMesh* ori_mesh;
	int m_size;
	std::vector<MeshProperty*> vProperties;
};