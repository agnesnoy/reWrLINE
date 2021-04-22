import numpy as np


# Maths


def cross(a, b):
    """
    Cross product of two arrays of 3-vectors
    """
    assert(np.shape(a)[0] == 3)
    assert(np.shape(a) == np.shape(b))
    return np.array([a[1]*b[2] - a[2]*b[1],
                     a[2]*b[0] - a[0]*b[2],
                     a[0]*b[1] - a[1]*b[0]])


def dot(a, b):
    """
    Dot product of two arrays of vectors
    """
    assert(np.shape(a) == np.shape(b))
    return sum(a[i] * b[i] for i in range(len(a)))


def norm(a):
    """
    Norms of an array of vectors
    """
    return np.sqrt(sum(x**2 for x in a))


def rotate_to_x(vector):
    """
    Returns the rotation matrix that rotates a vector to the x axis
    """
    normalised = vector / np.linalg.norm(vector)

    # C = -arctan(y / x)
    c = -np.arctan2(normalised[1], normalised[0])
    return np.array([[np.cos(c), -np.sin(c), 0.0],
                     [np.sin(c), np.cos(c),  0.0],
                     [0.0,       0.0,        1.0]])


def rotate_to_z(vector):
    """
    Returns the rotation matrix that rotates a vector to the z axis
    as follows:
      Define rotation matrix about (x, y) axis & Euler angles as
        Rxy(A, B) = Ry(B) `dot` Rx(A)
      Solve
        (0, 0, 1) = Rxy `dot` (x, y, z)
    """
    normalised = vector / np.linalg.norm(vector)

    # A = arctan(y / z)
    a = np.arctan2(normalised[1], normalised[2])
    # B = arctan(-x / |y, z|)
    # Explicitly NOT arctan2
    b = np.arctan(-normalised[0] /
                  np.sqrt(normalised[1]**2 + normalised[2]**2))

    rotate_x = np.array([[1.0,        0.0,       0.0],
                         [0.0,        np.cos(a), -np.sin(a)],
                         [0.0,        np.sin(a), np.cos(a)]])
    rotate_y = np.array([[np.cos(b),  0.0,       np.sin(b)],
                         [0.0,        1.0,       0.0],
                         [-np.sin(b), 0.0,       np.cos(b)]])

    return np.dot(rotate_y, rotate_x)


def twist(a1, a2, b1, b2, z):
    """
    Returns the twist
        a1 & a2 are the points defining the first vector
        b1 & b2 are the points defining the second vector
        z is the vector defining the z axis
    """
    vector_a = a2 - a1
    vector_b = b2 - b1

    unit_a = vector_a / np.linalg.norm(vector_a)
    unit_b = vector_b / np.linalg.norm(vector_b)
    unit_z = z / np.linalg.norm(z)

    # Generate rotation matrix rotating basis to the z axis
    rotation_matrix = rotate_to_z(unit_z)
    # Perform rotation
    unit_a = np.dot(rotation_matrix, unit_a)
    unit_b = np.dot(rotation_matrix, unit_b)

    # Now rotate such that a lies along the x axis
    rotated_b = np.dot(rotate_to_x(unit_a), unit_b)
    return np.arctan2(rotated_b[1], rotated_b[0]) * 180.0 / np.pi


# IO


def read(name, num_bp, num_steps, linear=False):
    """
    Reads a .mdcrd file & create returns a 3D array of atomic coordinates
    """
    with open(name + '/C.mdcrd', 'r') as file:
        file.readline()
        data = []
        for line in file:
            line_length = len(line)
            num_fields = line_length // 8
            for i in range(num_fields):
                data.append(float(line[i*8:(i+1)*8]))

        data_array = np.reshape(np.array(data), (len(data)//3, 3))

        x = np.reshape(data_array[:, 0], (num_steps, num_bp*2))
        y = np.reshape(data_array[:, 1], (num_steps, num_bp*2))
        z = np.reshape(data_array[:, 2], (num_steps, num_bp*2))

        strand_a = np.array([x[:, 0:num_bp], y[:, 0:num_bp], z[:, 0:num_bp]])
        strand_b = np.array([x[:, (2*num_bp - 1):(num_bp - 1):-1],
                             y[:, (2*num_bp - 1):(num_bp - 1):-1],
                             z[:, (2*num_bp - 1):(num_bp - 1):-1]])

        # Coordinate representation of a base-pair step
        midpoints = np.zeros(np.shape(strand_a))
        for i in range(num_bp):
            if linear and i+1 >= np.shape(strand_a)[2]:
                midpoints[:, :, i] = 0.5 * (strand_a[:, :, i] +
                                            strand_b[:, :, i])
            else:
                midpoints[:, :, i] = (0.25 *
                                      (strand_a[:, :, i] +
                                       strand_a[:, :, ((i + 1) % num_bp)] +
                                       strand_b[:, :, i] +
                                       strand_b[:, :, ((i + 1) % num_bp)]))

        return strand_a, strand_b, midpoints


def make_files(name, num_bp, num_steps, midpoints, caxis):
    """
    Makes xyz files of average C1' single helix
        midpoints is the array of midpoints of the C1' atoms of neighbouring
            base pairs
        caxis is the fully processed helical CAXIS
    """
    with open(name + '/C.xyz', 'w') as c_xyz, \
            open(name + '/C1.xyz', 'w') as c1_xyz, \
            open(name + '/C.3col', 'w') as c_3col, \
            open(name + '/C1.3col', 'w') as c1_3col:
        for i in range(num_steps):
            c_xyz.write(f"{num_bp}\n\n")
            c1_xyz.write(f"{num_bp}\n\n")
            print(f"\r\tStep {i}...", end=" ")
            for j in range(num_bp):
                c_xyz.write(" ".join(["H",
                                      f"{midpoints[0][i][j]:8.3f}",
                                      f"{midpoints[1][i][j]:8.3f}",
                                      f"{midpoints[2][i][j]:8.3f}",
                                      "\n"]))
                c1_xyz.write(" ".join(["H",
                                       f"{caxis[0][i][j]:8.3f}",
                                       f"{caxis[1][i][j]:8.3f}",
                                       f"{caxis[2][i][j]:8.3f}",
                                       "\n"]))
                c_3col.write(" ".join([f"{midpoints[0][i][j]:8.3f}",
                                       f"{midpoints[1][i][j]:8.3f}",
                                       f"{midpoints[2][i][j]:8.3f}",
                                       "\n"]))
                c1_3col.write(" ".join([f"{caxis[0][i][j]:8.3f}",
                                        f"{caxis[1][i][j]:8.3f}",
                                        f"{caxis[2][i][j]:8.3f}",
                                        "\n"]))
        print("Done!")


def sinreg(name, num_bp, num_steps, midpoints, caxis, write=True):
    """
    Calculates sine of register angles
    """
    result = np.zeros((num_steps, num_bp + 1))
    for j in range(num_bp):
        print(f"\r\tBase pair {j}...", end=" ")
        m = list(map(int, np.linspace(j-1, j+1, num=3) % num_bp))
        # Vectors bent on a plane
        v0 = caxis[:, :, m[1]] - caxis[:, :, m[0]]
        v1 = caxis[:, :, m[2]] - caxis[:, :, m[1]]
        plane_vector = cross(v1, v0)
        minor_groove = midpoints[:, :, m[0]] - caxis[:, :, m[0]]
        # sinreg[:, j+1] = |M <cross> C| / (|M| |C|)
        # where M = minor_groove & C = plane_vector
        result[:, j+1] = (norm(cross(minor_groove, plane_vector)) /
                          (norm(minor_groove) * norm(plane_vector)))
        # To obtain the sign of the result,
        # where positive is a minor groove pointing into the circle,
        # take the dot product of the minor groove with the unit normal vector
        for i in range(num_steps):
            if dot(v1 - v0, minor_groove)[i] < 0:
                result[i, j+1] = -result[i, j+1]
    result[:, 0] = np.linspace(0.01, 0.01*num_steps, num=num_steps)
    if write:
        np.savetxt(name + '/sinreg.ser', result, fmt='%8.3f')
    print("Done!")
    return result


def helix_axis(num_bp, num_steps, midpoints, strand_a, linear=False):
    """
    Calculates the first-order helical axis
    without taking the weight.
    Used for twist calculation.
    """
    result = np.zeros(np.shape(strand_a))
    for j in range(num_bp):
        print(f"\r\tBase pair {j}...", end=" ")
        # Summation of coordinates
        summation = np.zeros((3, num_steps))
        for t in range(num_steps):
            # Sum single helix position
            summation[:, t] += midpoints[:, t, j]
            k = 0
            while k < 5:
                k += 1
                if linear:
                    try:
                        summation[:, t] += (midpoints[:, t, j-k] +
                                            midpoints[:, t, j+k])
                    except IndexError:
                        k -= 1
                        break
                else:
                    summation[:, t] += (midpoints[:, t, (j-k) % num_bp] +
                                        midpoints[:, t, (j+k) % num_bp])
            # Average helix (almost full turn)
            result[:, t, j] = summation[:, t] / (2*k + 1)
    print("Done!")
    return result


def full_twist(name, num_bp, num_steps, strand_a, strand_b, haxis,
               linear=False, write=True):
    """
    Calculates twist
    """
    result = np.zeros((num_steps, num_bp))
    for j in range(num_bp):
        print(f"\r\tBase pair {j}...", end=" ")
        for t in range(num_steps):
            # Linear special cases
            if linear:
                # Does haxis[:, :, j+1] exist?
                if np.shape(haxis)[2] > j + 1:
                    # If haxis[:, :, j-1] doesn't exist,
                    # approximate using half the range
                    if j > 0:
                        z = haxis[:, t, j+1] - haxis[:, t, j-1]
                    else:
                        z = 2 * (haxis[:, t, j+1] - haxis[:, t, j])
                    result[t, j] = twist(strand_a[:, t, j],
                                         strand_b[:, t, j],
                                         strand_a[:, t, (j+1)],
                                         strand_b[:, t, (j+1)],
                                         z)
                # If not, just return zero
                # This may not be the best way to handle this
                else:
                    result[t, j] = 0
            else:
                z = haxis[:, t, (j+1) % num_bp] - haxis[:, t, (j-1) % num_bp]
                result[t, j] = twist(strand_a[:, t, j],
                                     strand_b[:, t, j],
                                     strand_a[:, t, (j+1) % num_bp],
                                     strand_b[:, t, (j+1) % num_bp],
                                     z)
    if write:
        np.savetxt(name + '/tw.ser', result, fmt='%8.3f')
    print("Done!")
    return result


def caxis(name, num_bp, num_steps, midpoints, tw, linear=False):
    """
    Calculates the central helical axis
    by performing the running average of each bp with its 2*k neighbours
    & including the weight of the excess base pair
    """
    result = np.zeros(np.shape(midpoints))
    for j in range(num_bp):
        print(f"\r\tBase pair {j}...", end=" ")
        total_twist = np.zeros(num_steps)
        summation = np.zeros((3, num_steps))
        for t in range(num_steps):
            total_twist[t] += tw[t, j]
            summation[:, t] += midpoints[:, t, j]
            k = 0
            # Find the point where two more flanking steps would make twist
            # exceed 360 degrees
            while total_twist[t] < 360.0:
                k += 1
                prev = total_twist[t]
                # In linear DNA, only go as far as the ends
                # This might not be the best approach
                if linear and (j-k < 0 or j+k >= num_bp):
                    break
                else:
                    total_twist[t] += (tw[t, (j-k) % num_bp] +
                                       tw[t, (j+k) % num_bp])
                    # Sum single helix position
                    summation[:, t] += (midpoints[:, t, (j-k) % num_bp] +
                                        midpoints[:, t, (j+k) % num_bp])
            # If the linear condition was met,
            # the twist must be less than 360 degrees,
            # so there's no need to remove the last two flanking steps
            if linear and (j-k < 0 or j+k >= num_bp):
                weight = 0
            else:
                # Add the flanks with weight < 1
                weight = (360.0 - prev) / (total_twist[t] - prev)
                summation[:, t] -= ((1 - weight) *
                                    (midpoints[:, t, (j-k) % num_bp] +
                                     midpoints[:, t, (j+k) % num_bp]))
            weight_3d = np.array([weight, weight, weight])
            result[:, t, j] = summation[:, t] / (2*(k + weight_3d) - 1)
    print("Done!")
    return result
