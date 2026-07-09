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
package org.eclipse.elk.core.debug.grandom.gRandom;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EEnum;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.emf.ecore.EReference;

/**
 * <!-- begin-user-doc -->
 * The <b>Package</b> for the model.
 * It contains accessors for the meta objects to represent
 * <ul>
 *   <li>each class,</li>
 *   <li>each feature of each class,</li>
 *   <li>each enum,</li>
 *   <li>and each data type</li>
 * </ul>
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomFactory
 * @model kind="package"
 * @generated
 */
public interface GRandomPackage extends EPackage
{
  /**
   * The package name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNAME = "gRandom";

  /**
   * The package namespace URI.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNS_URI = "http://www.eclipse.org/elk/core/debug/grandom/GRandom";

  /**
   * The package namespace name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNS_PREFIX = "gRandom";

  /**
   * The singleton instance of the package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  GRandomPackage eINSTANCE = org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl.init();

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.RandGraphImpl <em>Rand Graph</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.RandGraphImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getRandGraph()
   * @generated
   */
  int RAND_GRAPH = 0;

  /**
   * The feature id for the '<em><b>Configs</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int RAND_GRAPH__CONFIGS = 0;

  /**
   * The number of structural features of the '<em>Rand Graph</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int RAND_GRAPH_FEATURE_COUNT = 1;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.ConfigurationImpl <em>Configuration</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.ConfigurationImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getConfiguration()
   * @generated
   */
  int CONFIGURATION = 1;

  /**
   * The feature id for the '<em><b>Samples</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__SAMPLES = 0;

  /**
   * The feature id for the '<em><b>Form</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__FORM = 1;

  /**
   * The feature id for the '<em><b>Nodes</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__NODES = 2;

  /**
   * The feature id for the '<em><b>Edges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__EDGES = 3;

  /**
   * The feature id for the '<em><b>MW</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__MW = 4;

  /**
   * The feature id for the '<em><b>Max Width</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__MAX_WIDTH = 5;

  /**
   * The feature id for the '<em><b>MD</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__MD = 6;

  /**
   * The feature id for the '<em><b>Max Degree</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__MAX_DEGREE = 7;

  /**
   * The feature id for the '<em><b>PF</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__PF = 8;

  /**
   * The feature id for the '<em><b>Fraction</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__FRACTION = 9;

  /**
   * The feature id for the '<em><b>Hierarchy</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__HIERARCHY = 10;

  /**
   * The feature id for the '<em><b>Seed</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__SEED = 11;

  /**
   * The feature id for the '<em><b>Format</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__FORMAT = 12;

  /**
   * The feature id for the '<em><b>Filename</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION__FILENAME = 13;

  /**
   * The number of structural features of the '<em>Configuration</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int CONFIGURATION_FEATURE_COUNT = 14;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl <em>Hierarchy</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getHierarchy()
   * @generated
   */
  int HIERARCHY = 2;

  /**
   * The feature id for the '<em><b>Levels</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int HIERARCHY__LEVELS = 0;

  /**
   * The feature id for the '<em><b>Edges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int HIERARCHY__EDGES = 1;

  /**
   * The feature id for the '<em><b>Num Hierarch Nodes</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int HIERARCHY__NUM_HIERARCH_NODES = 2;

  /**
   * The feature id for the '<em><b>Cross Hierarch Rel</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int HIERARCHY__CROSS_HIERARCH_REL = 3;

  /**
   * The number of structural features of the '<em>Hierarchy</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int HIERARCHY_FEATURE_COUNT = 4;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl <em>Edges</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getEdges()
   * @generated
   */
  int EDGES = 3;

  /**
   * The feature id for the '<em><b>Density</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__DENSITY = 0;

  /**
   * The feature id for the '<em><b>Total</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__TOTAL = 1;

  /**
   * The feature id for the '<em><b>Relative</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__RELATIVE = 2;

  /**
   * The feature id for the '<em><b>Outbound</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__OUTBOUND = 3;

  /**
   * The feature id for the '<em><b>NEdges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__NEDGES = 4;

  /**
   * The feature id for the '<em><b>Labels</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__LABELS = 5;

  /**
   * The feature id for the '<em><b>Self Loops</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES__SELF_LOOPS = 6;

  /**
   * The number of structural features of the '<em>Edges</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int EDGES_FEATURE_COUNT = 7;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl <em>Nodes</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getNodes()
   * @generated
   */
  int NODES = 4;

  /**
   * The feature id for the '<em><b>NNodes</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES__NNODES = 0;

  /**
   * The feature id for the '<em><b>Ports</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES__PORTS = 1;

  /**
   * The feature id for the '<em><b>Labels</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES__LABELS = 2;

  /**
   * The feature id for the '<em><b>Size</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES__SIZE = 3;

  /**
   * The feature id for the '<em><b>Remove Isolated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES__REMOVE_ISOLATED = 4;

  /**
   * The number of structural features of the '<em>Nodes</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int NODES_FEATURE_COUNT = 5;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.SizeImpl <em>Size</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.SizeImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getSize()
   * @generated
   */
  int SIZE = 5;

  /**
   * The feature id for the '<em><b>Height</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int SIZE__HEIGHT = 0;

  /**
   * The feature id for the '<em><b>Width</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int SIZE__WIDTH = 1;

  /**
   * The number of structural features of the '<em>Size</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int SIZE_FEATURE_COUNT = 2;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl <em>Ports</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getPorts()
   * @generated
   */
  int PORTS = 6;

  /**
   * The feature id for the '<em><b>Labels</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS__LABELS = 0;

  /**
   * The feature id for the '<em><b>Re Use</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS__RE_USE = 1;

  /**
   * The feature id for the '<em><b>Size</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS__SIZE = 2;

  /**
   * The feature id for the '<em><b>Constraint</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS__CONSTRAINT = 3;

  /**
   * The feature id for the '<em><b>Flow</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS__FLOW = 4;

  /**
   * The number of structural features of the '<em>Ports</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int PORTS_FEATURE_COUNT = 5;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl <em>Flow</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFlow()
   * @generated
   */
  int FLOW = 7;

  /**
   * The feature id for the '<em><b>Flow Type</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int FLOW__FLOW_TYPE = 0;

  /**
   * The feature id for the '<em><b>Side</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int FLOW__SIDE = 1;

  /**
   * The feature id for the '<em><b>Amount</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int FLOW__AMOUNT = 2;

  /**
   * The number of structural features of the '<em>Flow</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int FLOW_FEATURE_COUNT = 3;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl <em>Double Quantity</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getDoubleQuantity()
   * @generated
   */
  int DOUBLE_QUANTITY = 8;

  /**
   * The feature id for the '<em><b>Quant</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__QUANT = 0;

  /**
   * The feature id for the '<em><b>Min</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__MIN = 1;

  /**
   * The feature id for the '<em><b>Min Max</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__MIN_MAX = 2;

  /**
   * The feature id for the '<em><b>Max</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__MAX = 3;

  /**
   * The feature id for the '<em><b>Mean</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__MEAN = 4;

  /**
   * The feature id for the '<em><b>Gaussian</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__GAUSSIAN = 5;

  /**
   * The feature id for the '<em><b>Stddv</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY__STDDV = 6;

  /**
   * The number of structural features of the '<em>Double Quantity</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int DOUBLE_QUANTITY_FEATURE_COUNT = 7;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Formats <em>Formats</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Formats
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFormats()
   * @generated
   */
  int FORMATS = 9;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Form <em>Form</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Form
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getForm()
   * @generated
   */
  int FORM = 10;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Side <em>Side</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Side
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getSide()
   * @generated
   */
  int SIDE = 11;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.FlowType <em>Flow Type</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.FlowType
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFlowType()
   * @generated
   */
  int FLOW_TYPE = 12;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType <em>Constraint Type</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType
   * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getConstraintType()
   * @generated
   */
  int CONSTRAINT_TYPE = 13;


  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.RandGraph <em>Rand Graph</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Rand Graph</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.RandGraph
   * @generated
   */
  EClass getRandGraph();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.debug.grandom.gRandom.RandGraph#getConfigs <em>Configs</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Configs</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.RandGraph#getConfigs()
   * @see #getRandGraph()
   * @generated
   */
  EReference getRandGraph_Configs();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration <em>Configuration</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Configuration</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration
   * @generated
   */
  EClass getConfiguration();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSamples <em>Samples</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Samples</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSamples()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_Samples();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getForm <em>Form</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Form</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getForm()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_Form();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getNodes <em>Nodes</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Nodes</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getNodes()
   * @see #getConfiguration()
   * @generated
   */
  EReference getConfiguration_Nodes();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getEdges <em>Edges</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Edges</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getEdges()
   * @see #getConfiguration()
   * @generated
   */
  EReference getConfiguration_Edges();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMW <em>MW</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>MW</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMW()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_MW();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxWidth <em>Max Width</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Max Width</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxWidth()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_MaxWidth();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMD <em>MD</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>MD</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMD()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_MD();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxDegree <em>Max Degree</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Max Degree</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxDegree()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_MaxDegree();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isPF <em>PF</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>PF</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isPF()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_PF();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFraction <em>Fraction</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Fraction</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFraction()
   * @see #getConfiguration()
   * @generated
   */
  EReference getConfiguration_Fraction();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getHierarchy <em>Hierarchy</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Hierarchy</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getHierarchy()
   * @see #getConfiguration()
   * @generated
   */
  EReference getConfiguration_Hierarchy();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSeed <em>Seed</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Seed</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSeed()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_Seed();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFormat <em>Format</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Format</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFormat()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_Format();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFilename <em>Filename</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Filename</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFilename()
   * @see #getConfiguration()
   * @generated
   */
  EAttribute getConfiguration_Filename();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy <em>Hierarchy</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Hierarchy</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy
   * @generated
   */
  EClass getHierarchy();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getLevels <em>Levels</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Levels</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getLevels()
   * @see #getHierarchy()
   * @generated
   */
  EReference getHierarchy_Levels();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getEdges <em>Edges</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Edges</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getEdges()
   * @see #getHierarchy()
   * @generated
   */
  EReference getHierarchy_Edges();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getNumHierarchNodes <em>Num Hierarch Nodes</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Num Hierarch Nodes</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getNumHierarchNodes()
   * @see #getHierarchy()
   * @generated
   */
  EReference getHierarchy_NumHierarchNodes();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getCrossHierarchRel <em>Cross Hierarch Rel</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Cross Hierarch Rel</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getCrossHierarchRel()
   * @see #getHierarchy()
   * @generated
   */
  EReference getHierarchy_CrossHierarchRel();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges <em>Edges</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Edges</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges
   * @generated
   */
  EClass getEdges();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isDensity <em>Density</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Density</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isDensity()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_Density();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isTotal <em>Total</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Total</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isTotal()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_Total();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isRelative <em>Relative</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Relative</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isRelative()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_Relative();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isOutbound <em>Outbound</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Outbound</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isOutbound()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_Outbound();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#getNEdges <em>NEdges</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>NEdges</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#getNEdges()
   * @see #getEdges()
   * @generated
   */
  EReference getEdges_NEdges();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isLabels <em>Labels</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Labels</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isLabels()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_Labels();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isSelfLoops <em>Self Loops</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Self Loops</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges#isSelfLoops()
   * @see #getEdges()
   * @generated
   */
  EAttribute getEdges_SelfLoops();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes <em>Nodes</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Nodes</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes
   * @generated
   */
  EClass getNodes();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getNNodes <em>NNodes</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>NNodes</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getNNodes()
   * @see #getNodes()
   * @generated
   */
  EReference getNodes_NNodes();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getPorts <em>Ports</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Ports</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getPorts()
   * @see #getNodes()
   * @generated
   */
  EReference getNodes_Ports();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes#isLabels <em>Labels</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Labels</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes#isLabels()
   * @see #getNodes()
   * @generated
   */
  EAttribute getNodes_Labels();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getSize <em>Size</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Size</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes#getSize()
   * @see #getNodes()
   * @generated
   */
  EReference getNodes_Size();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes#isRemoveIsolated <em>Remove Isolated</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Remove Isolated</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes#isRemoveIsolated()
   * @see #getNodes()
   * @generated
   */
  EAttribute getNodes_RemoveIsolated();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size <em>Size</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Size</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Size
   * @generated
   */
  EClass getSize();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getHeight <em>Height</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Height</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Size#getHeight()
   * @see #getSize()
   * @generated
   */
  EReference getSize_Height();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getWidth <em>Width</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Width</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Size#getWidth()
   * @see #getSize()
   * @generated
   */
  EReference getSize_Width();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports <em>Ports</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Ports</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports
   * @generated
   */
  EClass getPorts();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#isLabels <em>Labels</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Labels</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports#isLabels()
   * @see #getPorts()
   * @generated
   */
  EAttribute getPorts_Labels();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getReUse <em>Re Use</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Re Use</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports#getReUse()
   * @see #getPorts()
   * @generated
   */
  EReference getPorts_ReUse();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getSize <em>Size</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Size</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports#getSize()
   * @see #getPorts()
   * @generated
   */
  EReference getPorts_Size();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getConstraint <em>Constraint</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Constraint</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports#getConstraint()
   * @see #getPorts()
   * @generated
   */
  EAttribute getPorts_Constraint();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getFlow <em>Flow</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Flow</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports#getFlow()
   * @see #getPorts()
   * @generated
   */
  EReference getPorts_Flow();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow <em>Flow</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Flow</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Flow
   * @generated
   */
  EClass getFlow();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getFlowType <em>Flow Type</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Flow Type</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Flow#getFlowType()
   * @see #getFlow()
   * @generated
   */
  EAttribute getFlow_FlowType();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getSide <em>Side</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Side</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Flow#getSide()
   * @see #getFlow()
   * @generated
   */
  EAttribute getFlow_Side();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getAmount <em>Amount</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Amount</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Flow#getAmount()
   * @see #getFlow()
   * @generated
   */
  EReference getFlow_Amount();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity <em>Double Quantity</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Double Quantity</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity
   * @generated
   */
  EClass getDoubleQuantity();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getQuant <em>Quant</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Quant</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getQuant()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Quant();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMin <em>Min</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Min</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMin()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Min();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isMinMax <em>Min Max</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Min Max</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isMinMax()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_MinMax();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMax <em>Max</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Max</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMax()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Max();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMean <em>Mean</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Mean</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMean()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Mean();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isGaussian <em>Gaussian</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Gaussian</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isGaussian()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Gaussian();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getStddv <em>Stddv</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Stddv</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getStddv()
   * @see #getDoubleQuantity()
   * @generated
   */
  EAttribute getDoubleQuantity_Stddv();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.debug.grandom.gRandom.Formats <em>Formats</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Formats</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Formats
   * @generated
   */
  EEnum getFormats();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.debug.grandom.gRandom.Form <em>Form</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Form</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Form
   * @generated
   */
  EEnum getForm();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.debug.grandom.gRandom.Side <em>Side</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Side</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Side
   * @generated
   */
  EEnum getSide();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.debug.grandom.gRandom.FlowType <em>Flow Type</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Flow Type</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.FlowType
   * @generated
   */
  EEnum getFlowType();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType <em>Constraint Type</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Constraint Type</em>'.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType
   * @generated
   */
  EEnum getConstraintType();

  /**
   * Returns the factory that creates the instances of the model.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the factory that creates the instances of the model.
   * @generated
   */
  GRandomFactory getGRandomFactory();

  /**
   * <!-- begin-user-doc -->
   * Defines literals for the meta objects that represent
   * <ul>
   *   <li>each class,</li>
   *   <li>each feature of each class,</li>
   *   <li>each enum,</li>
   *   <li>and each data type</li>
   * </ul>
   * <!-- end-user-doc -->
   * @generated
   */
  interface Literals
  {
    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.RandGraphImpl <em>Rand Graph</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.RandGraphImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getRandGraph()
     * @generated
     */
    EClass RAND_GRAPH = eINSTANCE.getRandGraph();

    /**
     * The meta object literal for the '<em><b>Configs</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference RAND_GRAPH__CONFIGS = eINSTANCE.getRandGraph_Configs();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.ConfigurationImpl <em>Configuration</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.ConfigurationImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getConfiguration()
     * @generated
     */
    EClass CONFIGURATION = eINSTANCE.getConfiguration();

    /**
     * The meta object literal for the '<em><b>Samples</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__SAMPLES = eINSTANCE.getConfiguration_Samples();

    /**
     * The meta object literal for the '<em><b>Form</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__FORM = eINSTANCE.getConfiguration_Form();

    /**
     * The meta object literal for the '<em><b>Nodes</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference CONFIGURATION__NODES = eINSTANCE.getConfiguration_Nodes();

    /**
     * The meta object literal for the '<em><b>Edges</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference CONFIGURATION__EDGES = eINSTANCE.getConfiguration_Edges();

    /**
     * The meta object literal for the '<em><b>MW</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__MW = eINSTANCE.getConfiguration_MW();

    /**
     * The meta object literal for the '<em><b>Max Width</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__MAX_WIDTH = eINSTANCE.getConfiguration_MaxWidth();

    /**
     * The meta object literal for the '<em><b>MD</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__MD = eINSTANCE.getConfiguration_MD();

    /**
     * The meta object literal for the '<em><b>Max Degree</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__MAX_DEGREE = eINSTANCE.getConfiguration_MaxDegree();

    /**
     * The meta object literal for the '<em><b>PF</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__PF = eINSTANCE.getConfiguration_PF();

    /**
     * The meta object literal for the '<em><b>Fraction</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference CONFIGURATION__FRACTION = eINSTANCE.getConfiguration_Fraction();

    /**
     * The meta object literal for the '<em><b>Hierarchy</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference CONFIGURATION__HIERARCHY = eINSTANCE.getConfiguration_Hierarchy();

    /**
     * The meta object literal for the '<em><b>Seed</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__SEED = eINSTANCE.getConfiguration_Seed();

    /**
     * The meta object literal for the '<em><b>Format</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__FORMAT = eINSTANCE.getConfiguration_Format();

    /**
     * The meta object literal for the '<em><b>Filename</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute CONFIGURATION__FILENAME = eINSTANCE.getConfiguration_Filename();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl <em>Hierarchy</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getHierarchy()
     * @generated
     */
    EClass HIERARCHY = eINSTANCE.getHierarchy();

    /**
     * The meta object literal for the '<em><b>Levels</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference HIERARCHY__LEVELS = eINSTANCE.getHierarchy_Levels();

    /**
     * The meta object literal for the '<em><b>Edges</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference HIERARCHY__EDGES = eINSTANCE.getHierarchy_Edges();

    /**
     * The meta object literal for the '<em><b>Num Hierarch Nodes</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference HIERARCHY__NUM_HIERARCH_NODES = eINSTANCE.getHierarchy_NumHierarchNodes();

    /**
     * The meta object literal for the '<em><b>Cross Hierarch Rel</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference HIERARCHY__CROSS_HIERARCH_REL = eINSTANCE.getHierarchy_CrossHierarchRel();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl <em>Edges</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getEdges()
     * @generated
     */
    EClass EDGES = eINSTANCE.getEdges();

    /**
     * The meta object literal for the '<em><b>Density</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__DENSITY = eINSTANCE.getEdges_Density();

    /**
     * The meta object literal for the '<em><b>Total</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__TOTAL = eINSTANCE.getEdges_Total();

    /**
     * The meta object literal for the '<em><b>Relative</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__RELATIVE = eINSTANCE.getEdges_Relative();

    /**
     * The meta object literal for the '<em><b>Outbound</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__OUTBOUND = eINSTANCE.getEdges_Outbound();

    /**
     * The meta object literal for the '<em><b>NEdges</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference EDGES__NEDGES = eINSTANCE.getEdges_NEdges();

    /**
     * The meta object literal for the '<em><b>Labels</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__LABELS = eINSTANCE.getEdges_Labels();

    /**
     * The meta object literal for the '<em><b>Self Loops</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute EDGES__SELF_LOOPS = eINSTANCE.getEdges_SelfLoops();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl <em>Nodes</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getNodes()
     * @generated
     */
    EClass NODES = eINSTANCE.getNodes();

    /**
     * The meta object literal for the '<em><b>NNodes</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference NODES__NNODES = eINSTANCE.getNodes_NNodes();

    /**
     * The meta object literal for the '<em><b>Ports</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference NODES__PORTS = eINSTANCE.getNodes_Ports();

    /**
     * The meta object literal for the '<em><b>Labels</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute NODES__LABELS = eINSTANCE.getNodes_Labels();

    /**
     * The meta object literal for the '<em><b>Size</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference NODES__SIZE = eINSTANCE.getNodes_Size();

    /**
     * The meta object literal for the '<em><b>Remove Isolated</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute NODES__REMOVE_ISOLATED = eINSTANCE.getNodes_RemoveIsolated();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.SizeImpl <em>Size</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.SizeImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getSize()
     * @generated
     */
    EClass SIZE = eINSTANCE.getSize();

    /**
     * The meta object literal for the '<em><b>Height</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference SIZE__HEIGHT = eINSTANCE.getSize_Height();

    /**
     * The meta object literal for the '<em><b>Width</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference SIZE__WIDTH = eINSTANCE.getSize_Width();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl <em>Ports</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getPorts()
     * @generated
     */
    EClass PORTS = eINSTANCE.getPorts();

    /**
     * The meta object literal for the '<em><b>Labels</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute PORTS__LABELS = eINSTANCE.getPorts_Labels();

    /**
     * The meta object literal for the '<em><b>Re Use</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference PORTS__RE_USE = eINSTANCE.getPorts_ReUse();

    /**
     * The meta object literal for the '<em><b>Size</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference PORTS__SIZE = eINSTANCE.getPorts_Size();

    /**
     * The meta object literal for the '<em><b>Constraint</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute PORTS__CONSTRAINT = eINSTANCE.getPorts_Constraint();

    /**
     * The meta object literal for the '<em><b>Flow</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference PORTS__FLOW = eINSTANCE.getPorts_Flow();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl <em>Flow</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFlow()
     * @generated
     */
    EClass FLOW = eINSTANCE.getFlow();

    /**
     * The meta object literal for the '<em><b>Flow Type</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute FLOW__FLOW_TYPE = eINSTANCE.getFlow_FlowType();

    /**
     * The meta object literal for the '<em><b>Side</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute FLOW__SIDE = eINSTANCE.getFlow_Side();

    /**
     * The meta object literal for the '<em><b>Amount</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference FLOW__AMOUNT = eINSTANCE.getFlow_Amount();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl <em>Double Quantity</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getDoubleQuantity()
     * @generated
     */
    EClass DOUBLE_QUANTITY = eINSTANCE.getDoubleQuantity();

    /**
     * The meta object literal for the '<em><b>Quant</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__QUANT = eINSTANCE.getDoubleQuantity_Quant();

    /**
     * The meta object literal for the '<em><b>Min</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__MIN = eINSTANCE.getDoubleQuantity_Min();

    /**
     * The meta object literal for the '<em><b>Min Max</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__MIN_MAX = eINSTANCE.getDoubleQuantity_MinMax();

    /**
     * The meta object literal for the '<em><b>Max</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__MAX = eINSTANCE.getDoubleQuantity_Max();

    /**
     * The meta object literal for the '<em><b>Mean</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__MEAN = eINSTANCE.getDoubleQuantity_Mean();

    /**
     * The meta object literal for the '<em><b>Gaussian</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__GAUSSIAN = eINSTANCE.getDoubleQuantity_Gaussian();

    /**
     * The meta object literal for the '<em><b>Stddv</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute DOUBLE_QUANTITY__STDDV = eINSTANCE.getDoubleQuantity_Stddv();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Formats <em>Formats</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.Formats
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFormats()
     * @generated
     */
    EEnum FORMATS = eINSTANCE.getFormats();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Form <em>Form</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.Form
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getForm()
     * @generated
     */
    EEnum FORM = eINSTANCE.getForm();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Side <em>Side</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.Side
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getSide()
     * @generated
     */
    EEnum SIDE = eINSTANCE.getSide();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.FlowType <em>Flow Type</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.FlowType
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getFlowType()
     * @generated
     */
    EEnum FLOW_TYPE = eINSTANCE.getFlowType();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType <em>Constraint Type</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType
     * @see org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomPackageImpl#getConstraintType()
     * @generated
     */
    EEnum CONSTRAINT_TYPE = eINSTANCE.getConstraintType();

  }

} //GRandomPackage
