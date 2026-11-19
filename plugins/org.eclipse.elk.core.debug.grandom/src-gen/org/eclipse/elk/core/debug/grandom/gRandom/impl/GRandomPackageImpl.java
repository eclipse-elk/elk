/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom.impl;

import org.eclipse.elk.core.debug.grandom.gRandom.Configuration;
import org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType;
import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.Edges;
import org.eclipse.elk.core.debug.grandom.gRandom.Flow;
import org.eclipse.elk.core.debug.grandom.gRandom.FlowType;
import org.eclipse.elk.core.debug.grandom.gRandom.Form;
import org.eclipse.elk.core.debug.grandom.gRandom.Formats;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomFactory;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy;
import org.eclipse.elk.core.debug.grandom.gRandom.Nodes;
import org.eclipse.elk.core.debug.grandom.gRandom.Ports;
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph;
import org.eclipse.elk.core.debug.grandom.gRandom.Side;
import org.eclipse.elk.core.debug.grandom.gRandom.Size;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EEnum;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.emf.ecore.EReference;

import org.eclipse.emf.ecore.impl.EPackageImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model <b>Package</b>.
 * <!-- end-user-doc -->
 * @generated
 */
public class GRandomPackageImpl extends EPackageImpl implements GRandomPackage
{
  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass randGraphEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass configurationEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass hierarchyEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass edgesEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass nodesEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass sizeEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass portsEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass flowEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass doubleQuantityEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum formatsEEnum = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum formEEnum = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum sideEEnum = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum flowTypeEEnum = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum constraintTypeEEnum = null;

  /**
   * Creates an instance of the model <b>Package</b>, registered with
   * {@link org.eclipse.emf.ecore.EPackage.Registry EPackage.Registry} by the package
   * package URI value.
   * <p>Note: the correct way to create the package is via the static
   * factory method {@link #init init()}, which also performs
   * initialization of the package, or returns the registered package,
   * if one already exists.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.emf.ecore.EPackage.Registry
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#eNS_URI
   * @see #init()
   * @generated
   */
  private GRandomPackageImpl()
  {
    super(eNS_URI, GRandomFactory.eINSTANCE);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private static boolean isInited = false;

  /**
   * Creates, registers, and initializes the <b>Package</b> for this model, and for any others upon which it depends.
   *
   * <p>This method is used to initialize {@link GRandomPackage#eINSTANCE} when that field is accessed.
   * Clients should not invoke it directly. Instead, they should simply access that field to obtain the package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #eNS_URI
   * @see #createPackageContents()
   * @see #initializePackageContents()
   * @generated
   */
  public static GRandomPackage init()
  {
    if (isInited) return (GRandomPackage)EPackage.Registry.INSTANCE.getEPackage(GRandomPackage.eNS_URI);

    // Obtain or create and register package
    Object registeredGRandomPackage = EPackage.Registry.INSTANCE.get(eNS_URI);
    GRandomPackageImpl theGRandomPackage = registeredGRandomPackage instanceof GRandomPackageImpl ? (GRandomPackageImpl)registeredGRandomPackage : new GRandomPackageImpl();

    isInited = true;

    // Create package meta-data objects
    theGRandomPackage.createPackageContents();

    // Initialize created meta-data
    theGRandomPackage.initializePackageContents();

    // Mark meta-data to indicate it can't be changed
    theGRandomPackage.freeze();

    // Update the registry and return the package
    EPackage.Registry.INSTANCE.put(GRandomPackage.eNS_URI, theGRandomPackage);
    return theGRandomPackage;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getRandGraph()
  {
    return randGraphEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getRandGraph_Configs()
  {
    return (EReference)randGraphEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getConfiguration()
  {
    return configurationEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_Samples()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_Form()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getConfiguration_Nodes()
  {
    return (EReference)configurationEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getConfiguration_Edges()
  {
    return (EReference)configurationEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_MW()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_MaxWidth()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(5);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_MD()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(6);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_MaxDegree()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(7);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_PF()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(8);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getConfiguration_Fraction()
  {
    return (EReference)configurationEClass.getEStructuralFeatures().get(9);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getConfiguration_Hierarchy()
  {
    return (EReference)configurationEClass.getEStructuralFeatures().get(10);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_Seed()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(11);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_Format()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(12);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getConfiguration_Filename()
  {
    return (EAttribute)configurationEClass.getEStructuralFeatures().get(13);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getHierarchy()
  {
    return hierarchyEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getHierarchy_Levels()
  {
    return (EReference)hierarchyEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getHierarchy_Edges()
  {
    return (EReference)hierarchyEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getHierarchy_NumHierarchNodes()
  {
    return (EReference)hierarchyEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getHierarchy_CrossHierarchRel()
  {
    return (EReference)hierarchyEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getEdges()
  {
    return edgesEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_Density()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_Total()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_Relative()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_Outbound()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getEdges_NEdges()
  {
    return (EReference)edgesEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_Labels()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(5);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getEdges_SelfLoops()
  {
    return (EAttribute)edgesEClass.getEStructuralFeatures().get(6);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getNodes()
  {
    return nodesEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getNodes_NNodes()
  {
    return (EReference)nodesEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getNodes_Ports()
  {
    return (EReference)nodesEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getNodes_Labels()
  {
    return (EAttribute)nodesEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getNodes_Size()
  {
    return (EReference)nodesEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getNodes_RemoveIsolated()
  {
    return (EAttribute)nodesEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getSize()
  {
    return sizeEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getSize_Height()
  {
    return (EReference)sizeEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getSize_Width()
  {
    return (EReference)sizeEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getPorts()
  {
    return portsEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getPorts_Labels()
  {
    return (EAttribute)portsEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getPorts_ReUse()
  {
    return (EReference)portsEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getPorts_Size()
  {
    return (EReference)portsEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getPorts_Constraint()
  {
    return (EAttribute)portsEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getPorts_Flow()
  {
    return (EReference)portsEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getFlow()
  {
    return flowEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getFlow_FlowType()
  {
    return (EAttribute)flowEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getFlow_Side()
  {
    return (EAttribute)flowEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getFlow_Amount()
  {
    return (EReference)flowEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getDoubleQuantity()
  {
    return doubleQuantityEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Quant()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Min()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_MinMax()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Max()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Mean()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Gaussian()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(5);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getDoubleQuantity_Stddv()
  {
    return (EAttribute)doubleQuantityEClass.getEStructuralFeatures().get(6);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getFormats()
  {
    return formatsEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getForm()
  {
    return formEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getSide()
  {
    return sideEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getFlowType()
  {
    return flowTypeEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getConstraintType()
  {
    return constraintTypeEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public GRandomFactory getGRandomFactory()
  {
    return (GRandomFactory)getEFactoryInstance();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private boolean isCreated = false;

  /**
   * Creates the meta-model objects for the package.  This method is
   * guarded to have no affect on any invocation but its first.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void createPackageContents()
  {
    if (isCreated) return;
    isCreated = true;

    // Create classes and their features
    randGraphEClass = createEClass(RAND_GRAPH);
    createEReference(randGraphEClass, RAND_GRAPH__CONFIGS);

    configurationEClass = createEClass(CONFIGURATION);
    createEAttribute(configurationEClass, CONFIGURATION__SAMPLES);
    createEAttribute(configurationEClass, CONFIGURATION__FORM);
    createEReference(configurationEClass, CONFIGURATION__NODES);
    createEReference(configurationEClass, CONFIGURATION__EDGES);
    createEAttribute(configurationEClass, CONFIGURATION__MW);
    createEAttribute(configurationEClass, CONFIGURATION__MAX_WIDTH);
    createEAttribute(configurationEClass, CONFIGURATION__MD);
    createEAttribute(configurationEClass, CONFIGURATION__MAX_DEGREE);
    createEAttribute(configurationEClass, CONFIGURATION__PF);
    createEReference(configurationEClass, CONFIGURATION__FRACTION);
    createEReference(configurationEClass, CONFIGURATION__HIERARCHY);
    createEAttribute(configurationEClass, CONFIGURATION__SEED);
    createEAttribute(configurationEClass, CONFIGURATION__FORMAT);
    createEAttribute(configurationEClass, CONFIGURATION__FILENAME);

    hierarchyEClass = createEClass(HIERARCHY);
    createEReference(hierarchyEClass, HIERARCHY__LEVELS);
    createEReference(hierarchyEClass, HIERARCHY__EDGES);
    createEReference(hierarchyEClass, HIERARCHY__NUM_HIERARCH_NODES);
    createEReference(hierarchyEClass, HIERARCHY__CROSS_HIERARCH_REL);

    edgesEClass = createEClass(EDGES);
    createEAttribute(edgesEClass, EDGES__DENSITY);
    createEAttribute(edgesEClass, EDGES__TOTAL);
    createEAttribute(edgesEClass, EDGES__RELATIVE);
    createEAttribute(edgesEClass, EDGES__OUTBOUND);
    createEReference(edgesEClass, EDGES__NEDGES);
    createEAttribute(edgesEClass, EDGES__LABELS);
    createEAttribute(edgesEClass, EDGES__SELF_LOOPS);

    nodesEClass = createEClass(NODES);
    createEReference(nodesEClass, NODES__NNODES);
    createEReference(nodesEClass, NODES__PORTS);
    createEAttribute(nodesEClass, NODES__LABELS);
    createEReference(nodesEClass, NODES__SIZE);
    createEAttribute(nodesEClass, NODES__REMOVE_ISOLATED);

    sizeEClass = createEClass(SIZE);
    createEReference(sizeEClass, SIZE__HEIGHT);
    createEReference(sizeEClass, SIZE__WIDTH);

    portsEClass = createEClass(PORTS);
    createEAttribute(portsEClass, PORTS__LABELS);
    createEReference(portsEClass, PORTS__RE_USE);
    createEReference(portsEClass, PORTS__SIZE);
    createEAttribute(portsEClass, PORTS__CONSTRAINT);
    createEReference(portsEClass, PORTS__FLOW);

    flowEClass = createEClass(FLOW);
    createEAttribute(flowEClass, FLOW__FLOW_TYPE);
    createEAttribute(flowEClass, FLOW__SIDE);
    createEReference(flowEClass, FLOW__AMOUNT);

    doubleQuantityEClass = createEClass(DOUBLE_QUANTITY);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__QUANT);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__MIN);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__MIN_MAX);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__MAX);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__MEAN);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__GAUSSIAN);
    createEAttribute(doubleQuantityEClass, DOUBLE_QUANTITY__STDDV);

    // Create enums
    formatsEEnum = createEEnum(FORMATS);
    formEEnum = createEEnum(FORM);
    sideEEnum = createEEnum(SIDE);
    flowTypeEEnum = createEEnum(FLOW_TYPE);
    constraintTypeEEnum = createEEnum(CONSTRAINT_TYPE);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private boolean isInitialized = false;

  /**
   * Complete the initialization of the package and its meta-model.  This
   * method is guarded to have no affect on any invocation but its first.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void initializePackageContents()
  {
    if (isInitialized) return;
    isInitialized = true;

    // Initialize package
    setName(eNAME);
    setNsPrefix(eNS_PREFIX);
    setNsURI(eNS_URI);

    // Create type parameters

    // Set bounds for type parameters

    // Add supertypes to classes

    // Initialize classes and features; add operations and parameters
    initEClass(randGraphEClass, RandGraph.class, "RandGraph", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getRandGraph_Configs(), this.getConfiguration(), null, "configs", null, 0, -1, RandGraph.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(configurationEClass, Configuration.class, "Configuration", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getConfiguration_Samples(), ecorePackage.getEInt(), "samples", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_Form(), this.getForm(), "form", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getConfiguration_Nodes(), this.getNodes(), null, "nodes", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getConfiguration_Edges(), this.getEdges(), null, "edges", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_MW(), ecorePackage.getEBoolean(), "mW", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_MaxWidth(), ecorePackage.getEIntegerObject(), "maxWidth", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_MD(), ecorePackage.getEBoolean(), "mD", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_MaxDegree(), ecorePackage.getEIntegerObject(), "maxDegree", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_PF(), ecorePackage.getEBoolean(), "pF", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getConfiguration_Fraction(), this.getDoubleQuantity(), null, "fraction", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getConfiguration_Hierarchy(), this.getHierarchy(), null, "hierarchy", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_Seed(), ecorePackage.getEIntegerObject(), "seed", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_Format(), this.getFormats(), "format", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getConfiguration_Filename(), ecorePackage.getEString(), "filename", null, 0, 1, Configuration.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(hierarchyEClass, Hierarchy.class, "Hierarchy", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getHierarchy_Levels(), this.getDoubleQuantity(), null, "levels", null, 0, 1, Hierarchy.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getHierarchy_Edges(), this.getDoubleQuantity(), null, "edges", null, 0, 1, Hierarchy.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getHierarchy_NumHierarchNodes(), this.getDoubleQuantity(), null, "numHierarchNodes", null, 0, 1, Hierarchy.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getHierarchy_CrossHierarchRel(), this.getDoubleQuantity(), null, "crossHierarchRel", null, 0, 1, Hierarchy.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(edgesEClass, Edges.class, "Edges", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getEdges_Density(), ecorePackage.getEBoolean(), "density", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getEdges_Total(), ecorePackage.getEBoolean(), "total", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getEdges_Relative(), ecorePackage.getEBoolean(), "relative", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getEdges_Outbound(), ecorePackage.getEBoolean(), "outbound", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getEdges_NEdges(), this.getDoubleQuantity(), null, "nEdges", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getEdges_Labels(), ecorePackage.getEBoolean(), "labels", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getEdges_SelfLoops(), ecorePackage.getEBoolean(), "selfLoops", null, 0, 1, Edges.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(nodesEClass, Nodes.class, "Nodes", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getNodes_NNodes(), this.getDoubleQuantity(), null, "nNodes", null, 0, 1, Nodes.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getNodes_Ports(), this.getPorts(), null, "ports", null, 0, 1, Nodes.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getNodes_Labels(), ecorePackage.getEBoolean(), "labels", null, 0, 1, Nodes.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getNodes_Size(), this.getSize(), null, "size", null, 0, 1, Nodes.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getNodes_RemoveIsolated(), ecorePackage.getEBoolean(), "removeIsolated", null, 0, 1, Nodes.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(sizeEClass, Size.class, "Size", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getSize_Height(), this.getDoubleQuantity(), null, "height", null, 0, 1, Size.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getSize_Width(), this.getDoubleQuantity(), null, "width", null, 0, 1, Size.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(portsEClass, Ports.class, "Ports", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getPorts_Labels(), ecorePackage.getEBoolean(), "labels", null, 0, 1, Ports.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getPorts_ReUse(), this.getDoubleQuantity(), null, "reUse", null, 0, 1, Ports.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getPorts_Size(), this.getSize(), null, "size", null, 0, 1, Ports.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getPorts_Constraint(), this.getConstraintType(), "constraint", null, 0, 1, Ports.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getPorts_Flow(), this.getFlow(), null, "flow", null, 0, -1, Ports.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(flowEClass, Flow.class, "Flow", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getFlow_FlowType(), this.getFlowType(), "flowType", null, 0, 1, Flow.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getFlow_Side(), this.getSide(), "side", null, 0, 1, Flow.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getFlow_Amount(), this.getDoubleQuantity(), null, "amount", null, 0, 1, Flow.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(doubleQuantityEClass, DoubleQuantity.class, "DoubleQuantity", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getDoubleQuantity_Quant(), ecorePackage.getEDoubleObject(), "quant", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_Min(), ecorePackage.getEDoubleObject(), "min", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_MinMax(), ecorePackage.getEBoolean(), "minMax", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_Max(), ecorePackage.getEDoubleObject(), "max", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_Mean(), ecorePackage.getEDoubleObject(), "mean", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_Gaussian(), ecorePackage.getEBoolean(), "gaussian", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getDoubleQuantity_Stddv(), ecorePackage.getEDoubleObject(), "stddv", null, 0, 1, DoubleQuantity.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    // Initialize enums and add enum literals
    initEEnum(formatsEEnum, Formats.class, "Formats");
    addEEnumLiteral(formatsEEnum, Formats.ELKT);
    addEEnumLiteral(formatsEEnum, Formats.ELKG);

    initEEnum(formEEnum, Form.class, "Form");
    addEEnumLiteral(formEEnum, Form.TREES);
    addEEnumLiteral(formEEnum, Form.CUSTOM);
    addEEnumLiteral(formEEnum, Form.BIPARTITE);
    addEEnumLiteral(formEEnum, Form.BICONNECTED);
    addEEnumLiteral(formEEnum, Form.TRICONNECTED);
    addEEnumLiteral(formEEnum, Form.ACYCLIC);

    initEEnum(sideEEnum, Side.class, "Side");
    addEEnumLiteral(sideEEnum, Side.NORTH);
    addEEnumLiteral(sideEEnum, Side.EAST);
    addEEnumLiteral(sideEEnum, Side.SOUTH);
    addEEnumLiteral(sideEEnum, Side.WEST);

    initEEnum(flowTypeEEnum, FlowType.class, "FlowType");
    addEEnumLiteral(flowTypeEEnum, FlowType.INCOMING);
    addEEnumLiteral(flowTypeEEnum, FlowType.OUTGOING);

    initEEnum(constraintTypeEEnum, ConstraintType.class, "ConstraintType");
    addEEnumLiteral(constraintTypeEEnum, ConstraintType.FREE);
    addEEnumLiteral(constraintTypeEEnum, ConstraintType.SIDE);
    addEEnumLiteral(constraintTypeEEnum, ConstraintType.POSITION);
    addEEnumLiteral(constraintTypeEEnum, ConstraintType.ORDER);
    addEEnumLiteral(constraintTypeEEnum, ConstraintType.RATIO);

    // Create resource
    createResource(eNS_URI);
  }

} //GRandomPackageImpl
