#include "MeshProcessor.h"
#include <stdexcept>
#include <algorithm>
#include <cassert>


double& MeshFunction::operator[](int idx)
{
	return m_function.at(idx);
}

double MeshFunction::operator[] (int idx) const 
{
	return m_function[idx];
}

void MeshFunction::copyValues( const std::vector<double>& values )
{
	assert(m_size == values.size());
	m_function = values;
}

bool MeshProcessor::addProperty( MeshProperty* newProperty )
{
	vProperties.push_back(newProperty);
	return true;
}

MeshProcessor::MeshProcessor()
{
	mMesh = NULL;
}

void MeshProcessor::setMesh( CMesh* newMesh )
{
	mMesh = newMesh;
}

MeshProperty* MeshProcessor::retrievePropertyByID( int rid )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter) {
		if ((*iter)->id == rid)
			return *iter;
	}
	return NULL;
}

const MeshProperty* MeshProcessor::retrievePropertyByID( int rid ) const
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter) {
		if ((*iter)->id == rid)
			return *iter;
	}
	return NULL;
}

MeshProperty* MeshProcessor::retrievePropertyByName( const std::string& rn )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter) {
		if ((*iter)->name == rn)
			return *iter;
	}
	return NULL;
}

void MeshProcessor::removePropertyByID( int rid )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end();) {
		if ((*iter)->id == rid) {
			delete *iter;
			iter = vProperties.erase(iter);
		}
		else ++iter;
	}
}

void MeshProcessor::removePropertyByName( const std::string& rn )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter) {
		if ((*iter)->name == rn) {
			vProperties.erase(iter);
			break;
		}
	}
}

MeshProcessor::~MeshProcessor()
{
	for (MeshProperty* p : vProperties) delete p;	
}

void MeshProcessor::replaceProperty( MeshProperty* newProperty )
{
	removePropertyByID(newProperty->id);
	addProperty(newProperty);
}

MeshFeatureList::~MeshFeatureList()
{
	for (MeshFeature* f : m_vFeatures) delete f;
}
