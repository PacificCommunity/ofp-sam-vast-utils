#' Calculate the clockwise angle between two vectors
#' Vectors are specified relative to the origin
#' @param vector1 # vector of length 2. First entry is x and second entry is y
#' @param vector2 # vector of length 2. First entry is x and second entry is y
#' @return An angle in degrees
#' @export 

calc.clockwise.angle = function(vector1, vector2)
# https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors/16544330
{
	dot = vector1[1]*vector2[1] + vector1[2]*vector2[2]      # dot product between [x1, y1] and [x2, y2]
	det = vector1[1]*vector2[2] - vector1[2]*vector2[1]      # determinant
	angle = atan2(det, dot)*180/pi  # atan2(y, x) or atan2(sin, cos)
	if(angle<=0)
	{
		angle = abs(angle)
	} else {
		angle = 360 - angle
	}
	return(angle)
}
