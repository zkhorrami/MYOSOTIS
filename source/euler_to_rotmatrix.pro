Function euler_to_rotmatrix, euler_angles_in, DEGREE=degree

  ;Transform degrees to degree if asked for
  if keyword_set(degree) then $
    euler_angles = euler_angles_in * !pi/180. else $
    euler_angles = euler_angles_in

  ;Build three rotation matrices around the x-, y- and z-axis
  Rx = [[1.,0.,0.],$
    [0.,cos(euler_angles[0]),-sin(euler_angles[0])],$
    [0.,sin(euler_angles[0]), cos(euler_angles[0])]]
  Ry = [[cos(euler_angles[1]),0.,sin(euler_angles[1])],$
    [0,1.,0.],$
    [-sin(euler_angles[1]),0.,cos(euler_angles[1])]]
  Rz = [[cos(euler_angles[2]),-sin(euler_angles[2]),0.],$
    [sin(euler_angles[2]), cos(euler_angles[2]),0.],$
    [0.,0.,1.]]
  R = Rx#Ry#Rz

  ;Return the result
  return, R

End