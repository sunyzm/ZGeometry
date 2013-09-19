#include "ShapeEditor.h"


void ShapeEditor::copy( MeshCoordinates& coords ) const
{
	mProcessor->getMesh()->getVertCoordinates(coords);
}
