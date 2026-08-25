//==============================================================================
// FW_SquareChannel_CHT.geo
//
// Accepted first-wall geometry for CBS3D++_SI.
//
// Geometry:
//   streamwise length       = 200 mm
//   outer EUROFER section   = 21 mm x 21 mm
//   helium channel          = 15 mm x 10 mm
//   side wall thickness     = 3 mm
//   top wall thickness      = 3 mm
//   bottom wall thickness   = 8 mm
//   flow direction          = +x
//
// CBS physical IDs:
//   511 inlet
//   520 pressure outlet
//   530 external adiabatic solid wall
//   532 external TOP heat-flux wall
//    10 fluid volume
//    20 solid volume
//
// Mesh targets:
//   interface               = 0.15 mm
//   fluid interior          = 0.60 mm
//   solid bulk              = 1.00 mm
//   interface fine depth    = 0.20 mm
//   transition distance     = 1.50 mm
//
// Genuine 3-D HXT TET4 mesh. No Transfinite, Extrude, layers, prisms or
// prism-derived tetrahedral topology.
//==============================================================================

SetFactory("OpenCASCADE");
Geometry.OCCBoundsUseStl = 1;
Geometry.OCCBooleanPreserveNumbering = 1;

//==============================================================================
// 1. GEOMETRY PARAMETERS
//==============================================================================

L = 0.200;

Wout = 0.021;
Hout = 0.021;

Wfluid = 0.015;
Hfluid = 0.010;

tSide = 0.003;
tTop = 0.003;
tBottom = Hout - Hfluid - tTop;

yOuterMin = -Wout / 2.0;
yOuterMax =  Wout / 2.0;
zOuterMin = 0.0;
zOuterMax = Hout;

yFluidMin = -Wfluid / 2.0;
yFluidMax =  Wfluid / 2.0;
zFluidMin = tBottom;
zFluidMax = tBottom + Hfluid;

// OCC entity-selection tolerance only, not a mesh size.
eps = 1.0e-6;

//==============================================================================
// 2. CONFORMAL FLUID + SOLID VOLUMES
//==============================================================================

vOuter = newv;
Box(vOuter) =
{
    0.0, yOuterMin, zOuterMin,
    L,   Wout,      Hout
};

vFluid = newv;
Box(vFluid) =
{
    0.0, yFluidMin, zFluidMin,
    L,   Wfluid,    Hfluid
};

// Retain the fluid tool volume. Delete only the original outer block.
vSolid() = BooleanDifference
{
    Volume{vOuter};
    Delete;
}
{
    Volume{vFluid};
};

Coherence;

Physical Volume("fluid", 10) = {vFluid};
Physical Volume("solid", 20) = {vSolid()};

//==============================================================================
// 3. FOUR CONFORMAL FLUID-SOLID INTERFACE SURFACES
//    Used only for mesh sizing. Not exported as an external CBS boundary.
//==============================================================================

sInterfaceTop() = Surface In BoundingBox
{
    -eps, yFluidMin-eps, zFluidMax-eps,
    L+eps, yFluidMax+eps, zFluidMax+eps
};

sInterfaceBottom() = Surface In BoundingBox
{
    -eps, yFluidMin-eps, zFluidMin-eps,
    L+eps, yFluidMax+eps, zFluidMin+eps
};

sInterfaceLeft() = Surface In BoundingBox
{
    -eps, yFluidMin-eps, zFluidMin-eps,
    L+eps, yFluidMin+eps, zFluidMax+eps
};

sInterfaceRight() = Surface In BoundingBox
{
    -eps, yFluidMax-eps, zFluidMin-eps,
    L+eps, yFluidMax+eps, zFluidMax+eps
};

sInterface() = {};
sInterface() += sInterfaceTop();
sInterface() += sInterfaceBottom();
sInterface() += sInterfaceLeft();
sInterface() += sInterfaceRight();

If (#sInterface() != 4)
    Error("FW200: expected exactly 4 fluid-solid interface surfaces");
EndIf

//==============================================================================
// 4. FLUID INLET / OUTLET
//==============================================================================

sInlet() = Surface In BoundingBox
{
    -eps, yFluidMin-eps, zFluidMin-eps,
     eps, yFluidMax+eps, zFluidMax+eps
};

If (#sInlet() != 1)
    Error("FW200: expected exactly 1 fluid inlet surface");
EndIf

Physical Surface("inlet", 511) = {sInlet()};

sOutlet() = Surface In BoundingBox
{
    L-eps, yFluidMin-eps, zFluidMin-eps,
    L+eps, yFluidMax+eps, zFluidMax+eps
};

If (#sOutlet() != 1)
    Error("FW200: expected exactly 1 fluid outlet surface");
EndIf

Physical Surface("outlet", 520) = {sOutlet()};

//==============================================================================
// 5. EXTERNAL TOP HEAT-FLUX SURFACE
//==============================================================================

sHeatFlux() = Surface In BoundingBox
{
    -eps, yOuterMin-eps, zOuterMax-eps,
    L+eps, yOuterMax+eps, zOuterMax+eps
};

If (#sHeatFlux() != 1)
    Error("FW200: expected exactly 1 external top heat-flux surface");
EndIf

Physical Surface("heat_flux", 532) = {sHeatFlux()};

//==============================================================================
// 6. REMAINING EXTERNAL SOLID SURFACES = ADIABATIC
//==============================================================================

sOuterBottom() = Surface In BoundingBox
{
    -eps, yOuterMin-eps, zOuterMin-eps,
    L+eps, yOuterMax+eps, zOuterMin+eps
};

sOuterLeft() = Surface In BoundingBox
{
    -eps, yOuterMin-eps, zOuterMin-eps,
    L+eps, yOuterMin+eps, zOuterMax+eps
};

sOuterRight() = Surface In BoundingBox
{
    -eps, yOuterMax-eps, zOuterMin-eps,
    L+eps, yOuterMax+eps, zOuterMax+eps
};

sX0All() = Surface In BoundingBox
{
    -eps, yOuterMin-eps, zOuterMin-eps,
     eps, yOuterMax+eps, zOuterMax+eps
};

sXLAll() = Surface In BoundingBox
{
    L-eps, yOuterMin-eps, zOuterMin-eps,
    L+eps, yOuterMax+eps, zOuterMax+eps
};

sX0Solid() = sX0All();
sX0Solid() -= sInlet();

sXLSolid() = sXLAll();
sXLSolid() -= sOutlet();

sAdiabatic() = {};
sAdiabatic() += sOuterBottom();
sAdiabatic() += sOuterLeft();
sAdiabatic() += sOuterRight();
sAdiabatic() += sX0Solid();
sAdiabatic() += sXLSolid();

Physical Surface("adiabatic_wall", 530) = {sAdiabatic()};

//==============================================================================
// 7. MESH SIZE PARAMETERS
//==============================================================================

hInterface = 0.00015;
hBulk = 0.00100;
dFine = 0.00020;
dTransition = 0.00150;

// Interface distance field.
Field[1] = Distance;
Field[1].SurfacesList = {sInterface()};
Field[1].Sampling = 200;

// Interface refinement: 0.15 mm near interface -> 1.00 mm away.
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].SizeMin = hInterface;
Field[2].SizeMax = hBulk;
Field[2].DistMin = dFine;
Field[2].DistMax = dTransition;
Field[2].Sigmoid = 1;

// Fluid-core target/cap = 0.60 mm. This geometric box coincides with the fluid
// volume. The interface Distance field is finer and therefore wins near walls.
Field[3] = Box;
Field[3].VIn = 0.00060;
Field[3].VOut = hBulk;
Field[3].XMin = 0.0;
Field[3].XMax = L;
Field[3].YMin = yFluidMin;
Field[3].YMax = yFluidMax;
Field[3].ZMin = zFluidMin;
Field[3].ZMax = zFluidMax;
Field[3].Thickness = 0.00050;

Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;

// Do not let unrelated automatic size sources override the explicit fields.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeMin = hInterface;
Mesh.MeshSizeMax = hBulk;

//==============================================================================
// 8. TRUE 3-D LINEAR TETRAHEDRAL MESH
//==============================================================================

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 10;
Mesh.ElementOrder = 1;
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.SaveAll = 0;

//==============================================================================
// 9. AUDIT PRINTS
//==============================================================================

Printf("==============================================================");
Printf("FW200 TRUE 3-D TET CHT MESH");
Printf("L                  = %.6e m", L);
Printf("Outer WxH          = %.6e x %.6e m", Wout, Hout);
Printf("Fluid WxH          = %.6e x %.6e m", Wfluid, Hfluid);
Printf("Fluid z range      = %.6e .. %.6e m", zFluidMin, zFluidMax);
Printf("Interface h target = %.6e m", hInterface);
Printf("Fluid-core target  = %.6e m", 0.00060);
Printf("Solid bulk target  = %.6e m", hBulk);
Printf("Interface surfaces = %g", #sInterface());
Printf("BC511 inlet faces  = %g", #sInlet());
Printf("BC520 outlet faces = %g", #sOutlet());
Printf("BC532 heat faces   = %g", #sHeatFlux());
Printf("BC530 adiabatic    = %g", #sAdiabatic());
Printf("==============================================================");
