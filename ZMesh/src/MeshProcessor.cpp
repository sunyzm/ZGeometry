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

double MeshFunction::InnerProduct( const MeshFunction& f1, const MeshFunction& f2 )
{
	if (f1.m_size != f2.m_size)
		throw std::runtime_error("Inner product of incompatible manifold function");

	int dim = f1.m_size;
	double retval = 0.0;
	for (int i = 0; i < dim; ++i)
	{
		retval += f1.m_function[i] * f2.m_function[i];
	}

	return retval;
}

double MeshFunction::norm() const
{
	return MeshFunction::InnerProduct(*this, *this);
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
	mesh = NULL;
	m_size = 0;
}

void MeshProcessor::setMesh( CMesh* newMesh )
{
	mesh = newMesh;
	m_size = mesh->vertCount();
}

MeshProperty* MeshProcessor::retrievePropertyByID( int rid )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter)
	{
		if ((*iter)->id == rid)
			return *iter;
	}
	return NULL;
}

const MeshProperty* MeshProcessor::retrievePropertyByID_const( int rid ) const
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter)
	{
		if ((*iter)->id == rid)
			return *iter;
	}
	return NULL;
}

MeshProperty* MeshProcessor::retrievePropertyByName( const std::string& rn )
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter)
	{
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
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter)
	{
		if ((*iter)->name == rn)
		{
			vProperties.erase(iter);
			break;
		}
	}
}

MeshProcessor::~MeshProcessor()
{
	for (auto iter = vProperties.begin(); iter != vProperties.end(); ++iter)
	{
		delete *iter;
	}
}

void MeshProcessor::replaceProperty( MeshProperty* newProperty )
{
	removePropertyByID(newProperty->id);
	addProperty(newProperty);
}

MeshFeatureList::~MeshFeatureList()
{
	std::for_each(m_vFeatures.begin(), m_vFeatures.end(), [](MeshFeature* iter){ 
		delete iter;
	});
}
