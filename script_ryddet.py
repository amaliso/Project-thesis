from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
from mesh import*
executeOnCaeStartup()
#exec("nodes_lines.py")
Mdb()

import numpy as np
import math 


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

model_name = 'Model-1'
part_name = 'Part-1'
material_name = 'PA-CF'
E_mod = 4.3e3
pois = 0.25
dens = 1.17e-9
radius = 1

plates = 'Part-2'


s = mdb.models[model_name]
p = mdb.models[model_name].Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY) 



# Specify the input file name
input_file = 'lines.txt'

# Initialize an empty list to store the lines
all_lines = []

# Open the file for reading
with open(input_file, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace and split the line into two points
        line = line.strip()
        points = line.split(' - ')

        # Parse the two points as tuples and append them to the list
        point1_data, point2_data = points
        point1 = tuple(map(float, point1_data.strip('()').split(', ')))
        point2 = tuple(map(float, point2_data.strip('()').split(', ')))
        all_lines.append((point1, point2))



def create_lines(lines):
    for line in lines:
        p.WirePolyLine(points=(line), mergeType=IMPRINT, meshable=ON)
    #p.ReferencePoint(point = (20*math.cos(math.pi/8),0,85/2))  
    #v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    #p.ReferencePoint(point=v1[55])

create_lines(all_lines)





def define_material(mat_name, E_modulus, Poisson, density):
    s.Material(name=mat_name)
    s.materials[mat_name].Elastic(table=((E_modulus, Poisson), ))
    s.materials[mat_name].Density(table=((density, ), ))

define_material(mat_name=material_name, E_modulus=E_mod, Poisson=pois, density=dens)



def sets_by_box(zmax, zmin, set_name):
    tol=1e-5
    edges1 = p.edges.getByBoundingBox(zMin=zmin-tol, zMax=zmax+tol)
    defined_set = p.Set(edges=edges1, name=set_name)
    return defined_set, edges1

#sets_by_box(0,0,"bottom_set")
#top_set, top_edges =sets_by_box(50,50,"top_set") 
#whole,edges_whole = sets_by_box(50,0,"whole_set")

def set_by_cyl(center1, center2, radius, z_min, z_max, y_val):
    tol=1e-5
    all_edges = p.edges.getByBoundingCylinder(center1,center2, radius)
    whole_set = p.Set(edges = all_edges, name = "Whole_set")
    bottom_edges = p.edges.getByBoundingBox(zMin = z_min-tol, zMax=z_min+tol) 
    top_edges = p.edges.getByBoundingBox(zMin=z_max-tol, zMax=z_max+tol)
    
    upper_edges = p.vertices.getByBoundingBox(yMin=y_val-tol, yMax=y_val+tol)
    upper_set = p.Set(vertices = upper_edges, name = "Upper_set")
    
    lower_edges = p.vertices.getByBoundingBox(yMin= -y_val-tol, yMax= -y_val+tol)
    lower_set = p.Set(vertices = lower_edges, name = "Lower_set")
    #top_set = p.Set(edges = top_edges, name = "top_set")
    #bottom_set = p.Set(edges = bottom_edges, name = "bot_set")

    #middle_set = p.SetByBoolean(name=set_name, operation=DIFFERENCE, sets=(whole_set,top_set,bottom_set))
    return whole_set, upper_set, lower_set


c1 = (0,0,0)
c2 = (0,0,85)
r = 22
whole_set, upper_set, lower_set = set_by_cyl(c1,c2,r,0,85, 18.48)



def set_by_box(z_min, z_max, y_min, y_max):
    tol = 1e-5
    all_edges = p.edges.getByBoundingBox(zMin = z_min-tol, zMax=z_max+tol)
    


def assign_beam_ori(defined_set):
    region = defined_set
    p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(1.0,0.0,0.0))
    

assign_beam_ori(whole_set)
#assign_beam_ori(top_set)
#assign_beam_ori(bottom_set)

def section_assign(set, radius, material, profile_name, section_name):
    s.CircularProfile(name=profile_name, r=radius)
    s.BeamSection(name=section_name, integration=DURING_ANALYSIS, poissonRatio=0, profile=profile_name, material=material, temperatureVar=LINEAR, consistentMassMatrix=False)
    p.SectionAssignment(region=set, sectionName=section_name, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

section_assign(whole_set, radius=radius, material = material_name, profile_name='Circular', section_name='Section-1' )

def define_plates(part_name, sketch_p1, sketch_p2, ref_point):
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s1.rectangle(point1=(-30.0, 45.0), point2=(30.0, -45.0))
   # p1 = mdb.models[model_name].Part(name=part_name, dimensionality=THREE_D, type=DISCRETE_RIGID)
    #p1 = mdb.models[model_name].Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY) 
    p1 = mdb.models[model_name].Part(name=part_name, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p1.BaseShell(sketch=s1)
    p1.ReferencePoint(point=ref_point)

define_plates(plates, (-30,45), (30,-45), (0,0,-10))



def assembly():
    a1 = s.rootAssembly
    a1.DatumCsysByDefault(CARTESIAN)
    #a.Instance(name=instance_name, part=p, dependent=ON)
    
    
    p = mdb.models['Model-1'].parts['Part-1']
    a1.Instance(name='Part-1-1', part=p, dependent=ON)
    
    p = mdb.models['Model-1'].parts['Part-2']
    a1.Instance(name='Plate', part=p, dependent=ON)
    a1.translate(instanceList=('Plate', ), vector=(0.0, 18.48, 42.5))
    
    a1.DatumCsysByThreePoints(name='Datum csys-2', coordSysType=CARTESIAN, origin=(
    18.48,18.48,42.5), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
    
    a1 = mdb.models['Model-1'].rootAssembly
    a1.rotate(instanceList=('Plate', ), axisPoint=(0.0, 18.48, 42.5), 
    axisDirection=(10.0, 0.0, 0.0), angle=90.0)
    a1.RadialInstancePattern(instanceList=('Plate', ), 
    point=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), number=2, totalAngle=180.0) #Lage to plater

assembly()

def static_step(step_name):
    s.StaticStep(name=step_name, previous='Initial')
    mdb.models['Model-1'].steps['Step-1'].setValues(initialInc=0.1, maxInc=0.1)
    mdb.models['Model-1'].FieldOutputRequest(name='F-Output-2', createStepName='Step-1', variables=('U', 'SF'))

static_step(step_name='Step-1')



def interaction():
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
    0.3, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
    fraction=0.005, elasticSlipStiffness=None)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=OFF, 
    constraintEnforcementMethod=DEFAULT)

    
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['Plate'].faces
    side2Faces1 = s1.findAt(((10.0, 18.48, 57.5), ))
    region1=a.Surface(side2Faces=side2Faces1, name='m_Surf-3')
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.findAt(((7.65, 18.48, 70.0), ), ((7.65, 18.48, 50.0), ), 
    ((7.65, 18.48, 35.0), ), ((7.65, 18.48, 15.0), ), ((7.65, 18.48, 85.0), ),  
    ((7.65, 18.48, 0.0), ), ((-7.65, 18.48, 35.0), ), ((-7.65, 18.48, 15.0), ), 
    ((-7.65, 18.48, 0.0), ), ((-7.65, 18.48, 70.0), ), ((-7.65, 18.48, 50.0), ), 
    ((-7.65, 18.48, 85.0), ))
    region2=a.Set(vertices=verts1, name='s_Set-7')
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-3', 
    createStepName='Step-1', main=region1, secondary=region2, sliding=FINITE, 
    thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, 
    initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
    
    
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['Plate-rad-2'].faces
    side2Faces1 = s1.findAt(((-10.0, -18.48, 57.5), ))
    region1=a.Surface(side2Faces=side2Faces1, name='m_Surf-4')
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.findAt(((7.65, -18.48, 85.0), ), ((7.65, -18.48, 70.0), ),  
    ((-7.65, -18.48, 85.0), ), ((-7.65, -18.48, 70.0), ),  ((7.65, -18.48, 50.0), ), 
    ((7.65, -18.48, 35.0), ), ((-7.65, -18.48, 50.0), ), ((-7.65, -18.48, 35.0), ), 
    ((-7.65, -18.48, 15.0), ), ((7.65, -18.48, 15.0), ), ((7.65, -18.48, 0.0), ),  
    ((-7.65, -18.48, 0.0), ))
    region2=a.Set(vertices=verts1, name='s_Set-8')
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-4', 
    createStepName='Step-1', main=region1, secondary=region2, sliding=FINITE, 
    thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, 
    initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
    
    
interaction()





def loads(zmin,zmax,zload, set_name, instance_name, step_name,load_name):
    tol=1e-5
    a = s.rootAssembly
    e1 = a.instances[instance_name].edges
    edges1=e1.getByBoundingBox(zMin=zmin-tol, zMax=zmax+tol)
    region=a.Set(edges=edges1, name=set_name)
    s.LineLoad(name=load_name, createStepName=step_name, region=region, comp3=-zload)
    
#loads(60,60,1000,"load_set", "Instance-1",step_name="Static",load_name="Line load")

def boundary_cond(zmin,zmax,set_name, u, instance_name, bc_name, step_name):
    tol=1e-5
    a = s.rootAssembly
    e1 = a.instances[instance_name].edges
    edges1=e1.getByBoundingBox(zMin=zmin-tol, zMax=zmax+tol)
    region=a.Set(edges=edges1, name=set_name)
    s.DisplacementBC(name=bc_name, createStepName=step_name, region=region, 
    u1=u[0], u2=u[1], u3=u[2], ur1=u[3], ur2=u[4], ur3=u[5], amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

#u = [0,0,0,0,0,0]
#boundary_cond(0,0,"bound_set", u, instance_name='Instance-1', bc_name="Boundary-1", step_name='Static')

def boundarys():
    '''
    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['Plate-rad-2'].faces
    faces1 = f1.findAt(((-10.0, -18.48, 57.5), ))
    region = a.Set(faces=faces1, name='Set-1')
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Static', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['Plate'].referencePoints
    refPoints1=(r1[2], )
    region = a.Set(referencePoints=refPoints1, name='Set-4')
    mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='Static', 
    region=region, u1=0.0, u2=-10.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    '''
    
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['Plate'].referencePoints
    refPoints1=(r1[2], )
    region = a.Set(referencePoints=refPoints1, name='Set-4')
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Step-1', 
    region=region, u1=0.0, u2=-20.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['Plate-rad-2'].referencePoints
    refPoints1=(r1[2], )
    region = a.Set(referencePoints=refPoints1, name='Set-5')
    mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.findAt(((-7.65, -18.48, 85.0), ))
    region = a.Set(vertices=verts1, name='Set-6')
    mdb.models['Model-1'].DisplacementBC(name='BC-3', createStepName='Step-1', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.findAt(((7.65, -18.48, 0.0), ))
    region = a.Set(vertices=verts1, name='Set-7')
    mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='Step-1', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
    
boundarys()

def mesh():
    p = mdb.models['Model-1'].parts['Part-2']
    p.seedPart(size=10.0, deviationFactor=0.1, minSizeValue=0.9)
    p = mdb.models['Model-1'].parts['Part-2']
    p.generateMesh()
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.seedPart(size=1.5, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['Part-1']
    p.generateMesh()

mesh()
#mesh(5.0, B31, whole)


def create_job(job_name, model_name):
    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

    mdb.jobs[job_name].submit(consistencyChecking=OFF)

#create_job()