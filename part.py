from OCC.Core.GProp import GProp_PEquation
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.TopExp import topexp, TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods, TopoDS_Edge, TopoDS_Face
from OCC.Core.BRep import BRep_Tool
import numpy as np
from scipy.spatial import ConvexHull
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.GeomAbs import GeomAbs_Cylinder
import math
from OCC.Core.TopTools import (
    TopTools_IndexedMapOfShape,
    TopTools_IndexedDataMapOfShapeListOfShape,
    TopTools_ListIteratorOfListOfShape
)
from OCC.Core.GCPnts import GCPnts_AbscissaPoint
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib


class Part:
    # properties object
    props = GProp_GProps()

    def __init__(self, shape):
        self.shape = shape

    # Calculate and print bounding box dimensions
    def getBoundingBox(self):
        vertices = []
        exp = TopExp_Explorer(self.shape, TopAbs_VERTEX)
        while exp.More():
            vertex = exp.Current()
            pnt = BRep_Tool.Pnt(vertex)
            vertices.append(pnt)
            exp.Next()

        # Convert list of gp_Pnt to TColgp_Array1OfPnt
        points_array = TColgp_Array1OfPnt(1, len(vertices))
        for i, p in enumerate(vertices, 1):
            points_array.SetValue(i, p)

        # Calculate the bounding box
        precision = 1e-6
        prop_eq = GProp_PEquation(points_array, precision)

        P = gp_Pnt()
        V1 = gp_Vec()
        V2 = gp_Vec()
        V3 = gp_Vec()
        prop_eq.Box(P, V1, V2, V3)

        return P, V1, V2, V3

    def getMaxDimMM(self):
        P, V1, V2, V3 = self.getBoundingBox()
        return max(V1.Magnitude(), V2.Magnitude(), V3.Magnitude())

    # Calculate volume from brepgrprops
    def getVolumeMM3(self):
        tolerance = 1e-5
        brepgprop.VolumeProperties(self.shape, self.props, tolerance)

        volume = self.props.Mass()
        return volume

    # Calculate surface area from brepgprops
    def getSurfaceAreaMM2(self):
        brepgprop.SurfaceProperties(self.shape, self.props)

        surface_area = self.props.Mass()
        return surface_area

    # Mesh the shape to get vertices
    def getVertices(self, mesh_deflection=0.5):
        BRepMesh_IncrementalMesh(self.shape, mesh_deflection)
        verts = []
        exp = TopExp_Explorer(self.shape, TopAbs_VERTEX)
        while exp.More():
            v = topods.Vertex(exp.Current())
            pnt = BRep_Tool.Pnt(v)
            verts.append([pnt.X(), pnt.Y(), pnt.Z()])
            exp.Next()
        return np.array(verts)

    # Get the volume of the smallest convex hull around the shape
    def getConvexHullVolumeMM3(self):
        vertices = self.getVertices()
        hull = ConvexHull(vertices)
        convex_hull_volume = hull.volume
        return convex_hull_volume



    def getNumHoles(self, min_radius=0.5, max_radius=50.0):
        seen = []
        holes = 0
        exp = TopExp_Explorer(self.shape, TopAbs_FACE)
        while exp.More():
            face = exp.Current()
            surf = BRepAdaptor_Surface(face)
            if surf.GetType() == GeomAbs_Cylinder:
                cyl = surf.Cylinder()
                r = cyl.Radius()
                axis = cyl.Axis().Direction()
                if min_radius <= r <= max_radius:
                        # deduplicate by (r, axis location)
                        key = (round(r, 2), round(cyl.Location().Z(), 2))
                        if key not in seen:
                            holes += 1
                            seen.append(key)
            exp.Next()
        return holes

    def getNumSharpEdges(self, angle_tol_deg=90, debug=False):
        """
        Count sharp edges of a shape using dihedral angle.
        angle_tol_deg = threshold for sharpness (smaller = stricter).
        """
        sharp = 0
        seen = TopTools_IndexedMapOfShape()
        edge_face_map = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp.MapShapesAndAncestors(self.shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map)

        exp = TopExp_Explorer(self.shape, TopAbs_EDGE)
        while exp.More():
            edge = topods.Edge(exp.Current())

            # deduplicate
            if seen.Contains(edge):
                exp.Next()
                continue
            seen.Add(edge)

            # filter out microscopic edges
            diag = bbox_diag(self.shape)
            min_edge_len = diag * 0.1
            #min_edge_len = 2
            try:
                L = edge_length(edge)
            except:
                L = 0
            if L < min_edge_len:
                exp.Next()
                continue

            faces = []
            face_list = edge_face_map.FindFromKey(edge)
            if not face_list.Size() == 0:
                it = TopTools_ListIteratorOfListOfShape(face_list)
                while it.More():
                    faces.append(topods.Face(it.Value()))
                    it.Next()

            if len(faces) == 2:
                angle = angle_between_normals(faces[0], faces[1], edge)
                if angle is not None and angle < angle_tol_deg:
                    sharp += 1
                    if debug:
                        print(f"Sharp edge found, angle={angle:.1f}")

            exp.Next()

        return sharp

    def getNumPockets(self, tol=0.5):
        # Simplified: count cylindrical faces with normals pointing inward
        pockets = []
        exp = TopExp_Explorer(self.shape, TopAbs_FACE)

        while exp.More():
            face = exp.Current()
            surf = BRepAdaptor_Surface(face)
            stype = surf.GetType()

            if stype == GeomAbs_Cylinder:
                u_min, u_max = surf.FirstUParameter(), surf.LastUParameter()
                v_min, v_max = surf.FirstVParameter(), surf.LastVParameter()
                u_mid, v_mid = (u_min + u_max) / 2, (v_min + v_max) / 2

                normal = surface_normal(face, u_mid, v_mid)
                if normal is None:
                    exp.Next()
                    continue

                pnt = surf.Value(u_mid, v_mid)
                vec_to_center = gp_Dir(-pnt.X(), -pnt.Y(), -pnt.Z())

                if normal.Dot(vec_to_center) > 0:  # inward
                    # Check for duplicates
                    duplicate = False
                    for existing in pockets:
                        dist = gp_Vec(pnt, existing).Magnitude()
                        if dist < tol:
                            duplicate = True
                            break
                    if not duplicate:
                        pockets.append(pnt)

            exp.Next()

        return len(pockets)

    def getComplexity(self):
        surface_area = self.getSurfaceAreaMM2()
        volume = self.getVolumeMM3()
        diag = bbox_diag(self.shape)
        size_scale = max(diag, 1.0)  # prevent division by zero

        r_star = (surface_area / (volume ** (2.0 / 3.0))) if volume > 0 else 0.0
        r_star_scaled = r_star * 0.1  # scale down to prevent tiny parts from dominating

        H = self.getNumHoles()
        print(f"\nHoles: {H}")
        E = self.getNumSharpEdges()
        print(f"Sharp Edges: {E}")
        P = self.getNumPockets()
        print(f"Pockets: {P}")

        # Set alpha weights based on part size
        if diag <= 50:  # small part (tiny gears)
            alpha_H, alpha_E, alpha_P, alpha_r = 0.1, 0.07, 0.1, 0.01
        elif diag <= 100:  # medium part
            alpha_H, alpha_E, alpha_P, alpha_r = 0.2, 0.1, 0.2, 0.02
        else:  # large part
            alpha_H, alpha_E, alpha_P, alpha_r = 1.5, 0.75, 1.5, 0.25

        # Weighted sum
        score_raw = (alpha_H * H / size_scale +
                     alpha_E * E / size_scale +
                     alpha_P * P / size_scale +
                     alpha_r * r_star_scaled)

        # Squash to 0–1
        complexity = math.tanh(score_raw)
        return complexity

def bbox_diag(shape):
    box = Bnd_Box()
    brepbndlib.Add(shape, box)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return ((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2) ** 0.5

def edge_length(edge):
    curve = BRepAdaptor_Curve(edge)
    return GCPnts_AbscissaPoint.Length(curve)

#Return a gp_Dir normal at param (u,v) on a face.
def surface_normal(face, u, v):
    surf = BRepAdaptor_Surface(face)

    P = gp_Pnt()
    D1U = gp_Vec()
    D1V = gp_Vec()

    surf.D1(u, v, P, D1U, D1V)   # fills P, D1U, D1V

    n = D1U.Crossed(D1V)
    if n.Magnitude() == 0:
        return None
    return gp_Dir(n)

def angle_between_normals(face1, face2, edge):
    # Curve on face1
    c2d_1, first1, last1 = BRep_Tool.CurveOnSurface(edge, face1)
    if c2d_1 is None:
        return None
    midparam1 = 0.5 * (first1 + last1)
    uv1 = c2d_1.Value(midparam1)

    # Curve on face2
    c2d_2, first2, last2 = BRep_Tool.CurveOnSurface(edge, face2)
    if c2d_2 is None:
        return None
    midparam2 = 0.5 * (first2 + last2)
    uv2 = c2d_2.Value(midparam2)

    # Compute normals
    n1 = surface_normal(face1, uv1.X(), uv1.Y())
    n2 = surface_normal(face2, uv2.X(), uv2.Y())
    if not n1 or not n2:
        return None

    # Angle between normals in degrees
    angle_rad = n1.Angle(n2)
    return math.degrees(angle_rad)