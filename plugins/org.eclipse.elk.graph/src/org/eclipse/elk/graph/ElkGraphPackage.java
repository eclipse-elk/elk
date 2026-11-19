/**
 * Copyright (c) 2016 Kiel University and others.
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
 */
package org.eclipse.elk.graph;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EDataType;
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
 * @see org.eclipse.elk.graph.ElkGraphFactory
 * @model kind="package"
 * @generated
 */
public interface ElkGraphPackage extends EPackage {
    /**
     * The package name.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    String eNAME = "graph";

    /**
     * The package namespace URI.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    String eNS_URI = "http://www.eclipse.org/elk/ElkGraph";

    /**
     * The package namespace name.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    String eNS_PREFIX = "graph";

    /**
     * The singleton instance of the package.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    ElkGraphPackage eINSTANCE = org.eclipse.elk.graph.impl.ElkGraphPackageImpl.init();

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.properties.IPropertyHolder <em>IProperty Holder</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.properties.IPropertyHolder
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getIPropertyHolder()
     * @generated
     */
    int IPROPERTY_HOLDER = 0;

    /**
     * The number of structural features of the '<em>IProperty Holder</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int IPROPERTY_HOLDER_FEATURE_COUNT = 0;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.EMapPropertyHolderImpl <em>EMap Property Holder</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.EMapPropertyHolderImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getEMapPropertyHolder()
     * @generated
     */
    int EMAP_PROPERTY_HOLDER = 1;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int EMAP_PROPERTY_HOLDER__PROPERTIES = IPROPERTY_HOLDER_FEATURE_COUNT + 0;

    /**
     * The number of structural features of the '<em>EMap Property Holder</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int EMAP_PROPERTY_HOLDER_FEATURE_COUNT = IPROPERTY_HOLDER_FEATURE_COUNT + 1;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkGraphElementImpl <em>Element</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkGraphElementImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkGraphElement()
     * @generated
     */
    int ELK_GRAPH_ELEMENT = 2;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_GRAPH_ELEMENT__PROPERTIES = EMAP_PROPERTY_HOLDER__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_GRAPH_ELEMENT__LABELS = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_GRAPH_ELEMENT__IDENTIFIER = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 1;

    /**
     * The number of structural features of the '<em>Element</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_GRAPH_ELEMENT_FEATURE_COUNT = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 2;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkShapeImpl <em>Elk Shape</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkShapeImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkShape()
     * @generated
     */
    int ELK_SHAPE = 3;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__PROPERTIES = ELK_GRAPH_ELEMENT__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__LABELS = ELK_GRAPH_ELEMENT__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__IDENTIFIER = ELK_GRAPH_ELEMENT__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Height</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__HEIGHT = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Width</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__WIDTH = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 1;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__X = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 2;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE__Y = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 3;

    /**
     * The number of structural features of the '<em>Elk Shape</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_SHAPE_FEATURE_COUNT = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 4;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkLabelImpl <em>Elk Label</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkLabelImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkLabel()
     * @generated
     */
    int ELK_LABEL = 4;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__PROPERTIES = ELK_SHAPE__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__LABELS = ELK_SHAPE__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__IDENTIFIER = ELK_SHAPE__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Height</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__HEIGHT = ELK_SHAPE__HEIGHT;

    /**
     * The feature id for the '<em><b>Width</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__WIDTH = ELK_SHAPE__WIDTH;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__X = ELK_SHAPE__X;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__Y = ELK_SHAPE__Y;

    /**
     * The feature id for the '<em><b>Parent</b></em>' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__PARENT = ELK_SHAPE_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Text</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL__TEXT = ELK_SHAPE_FEATURE_COUNT + 1;

    /**
     * The number of structural features of the '<em>Elk Label</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_LABEL_FEATURE_COUNT = ELK_SHAPE_FEATURE_COUNT + 2;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkConnectableShapeImpl <em>Elk Connectable Shape</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkConnectableShapeImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkConnectableShape()
     * @generated
     */
    int ELK_CONNECTABLE_SHAPE = 5;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__PROPERTIES = ELK_SHAPE__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__LABELS = ELK_SHAPE__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__IDENTIFIER = ELK_SHAPE__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Height</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__HEIGHT = ELK_SHAPE__HEIGHT;

    /**
     * The feature id for the '<em><b>Width</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__WIDTH = ELK_SHAPE__WIDTH;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__X = ELK_SHAPE__X;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__Y = ELK_SHAPE__Y;

    /**
     * The feature id for the '<em><b>Outgoing Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES = ELK_SHAPE_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Incoming Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE__INCOMING_EDGES = ELK_SHAPE_FEATURE_COUNT + 1;

    /**
     * The number of structural features of the '<em>Elk Connectable Shape</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_CONNECTABLE_SHAPE_FEATURE_COUNT = ELK_SHAPE_FEATURE_COUNT + 2;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkNodeImpl <em>Elk Node</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkNodeImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkNode()
     * @generated
     */
    int ELK_NODE = 6;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__PROPERTIES = ELK_CONNECTABLE_SHAPE__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__LABELS = ELK_CONNECTABLE_SHAPE__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__IDENTIFIER = ELK_CONNECTABLE_SHAPE__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Height</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__HEIGHT = ELK_CONNECTABLE_SHAPE__HEIGHT;

    /**
     * The feature id for the '<em><b>Width</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__WIDTH = ELK_CONNECTABLE_SHAPE__WIDTH;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__X = ELK_CONNECTABLE_SHAPE__X;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__Y = ELK_CONNECTABLE_SHAPE__Y;

    /**
     * The feature id for the '<em><b>Outgoing Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__OUTGOING_EDGES = ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES;

    /**
     * The feature id for the '<em><b>Incoming Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__INCOMING_EDGES = ELK_CONNECTABLE_SHAPE__INCOMING_EDGES;

    /**
     * The feature id for the '<em><b>Ports</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__PORTS = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Children</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__CHILDREN = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 1;

    /**
     * The feature id for the '<em><b>Parent</b></em>' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__PARENT = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 2;

    /**
     * The feature id for the '<em><b>Contained Edges</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__CONTAINED_EDGES = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 3;

    /**
     * The feature id for the '<em><b>Hierarchical</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE__HIERARCHICAL = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 4;

    /**
     * The number of structural features of the '<em>Elk Node</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_NODE_FEATURE_COUNT = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 5;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkPortImpl <em>Elk Port</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkPortImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkPort()
     * @generated
     */
    int ELK_PORT = 7;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__PROPERTIES = ELK_CONNECTABLE_SHAPE__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__LABELS = ELK_CONNECTABLE_SHAPE__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__IDENTIFIER = ELK_CONNECTABLE_SHAPE__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Height</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__HEIGHT = ELK_CONNECTABLE_SHAPE__HEIGHT;

    /**
     * The feature id for the '<em><b>Width</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__WIDTH = ELK_CONNECTABLE_SHAPE__WIDTH;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__X = ELK_CONNECTABLE_SHAPE__X;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__Y = ELK_CONNECTABLE_SHAPE__Y;

    /**
     * The feature id for the '<em><b>Outgoing Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__OUTGOING_EDGES = ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES;

    /**
     * The feature id for the '<em><b>Incoming Edges</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__INCOMING_EDGES = ELK_CONNECTABLE_SHAPE__INCOMING_EDGES;

    /**
     * The feature id for the '<em><b>Parent</b></em>' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT__PARENT = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 0;

    /**
     * The number of structural features of the '<em>Elk Port</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PORT_FEATURE_COUNT = ELK_CONNECTABLE_SHAPE_FEATURE_COUNT + 1;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkEdgeImpl <em>Elk Edge</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkEdgeImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkEdge()
     * @generated
     */
    int ELK_EDGE = 8;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__PROPERTIES = ELK_GRAPH_ELEMENT__PROPERTIES;

    /**
     * The feature id for the '<em><b>Labels</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__LABELS = ELK_GRAPH_ELEMENT__LABELS;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__IDENTIFIER = ELK_GRAPH_ELEMENT__IDENTIFIER;

    /**
     * The feature id for the '<em><b>Containing Node</b></em>' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__CONTAINING_NODE = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Sources</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__SOURCES = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 1;

    /**
     * The feature id for the '<em><b>Targets</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__TARGETS = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 2;

    /**
     * The feature id for the '<em><b>Sections</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__SECTIONS = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 3;

    /**
     * The feature id for the '<em><b>Hyperedge</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__HYPEREDGE = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 4;

    /**
     * The feature id for the '<em><b>Hierarchical</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__HIERARCHICAL = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 5;

    /**
     * The feature id for the '<em><b>Selfloop</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__SELFLOOP = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 6;

    /**
     * The feature id for the '<em><b>Connected</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE__CONNECTED = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 7;

    /**
     * The number of structural features of the '<em>Elk Edge</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_FEATURE_COUNT = ELK_GRAPH_ELEMENT_FEATURE_COUNT + 8;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkBendPointImpl <em>Elk Bend Point</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkBendPointImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkBendPoint()
     * @generated
     */
    int ELK_BEND_POINT = 9;

    /**
     * The feature id for the '<em><b>X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_BEND_POINT__X = 0;

    /**
     * The feature id for the '<em><b>Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_BEND_POINT__Y = 1;

    /**
     * The number of structural features of the '<em>Elk Bend Point</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_BEND_POINT_FEATURE_COUNT = 2;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl <em>Elk Edge Section</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkEdgeSectionImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkEdgeSection()
     * @generated
     */
    int ELK_EDGE_SECTION = 10;

    /**
     * The feature id for the '<em><b>Properties</b></em>' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__PROPERTIES = EMAP_PROPERTY_HOLDER__PROPERTIES;

    /**
     * The feature id for the '<em><b>Start X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__START_X = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 0;

    /**
     * The feature id for the '<em><b>Start Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__START_Y = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 1;

    /**
     * The feature id for the '<em><b>End X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__END_X = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 2;

    /**
     * The feature id for the '<em><b>End Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__END_Y = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 3;

    /**
     * The feature id for the '<em><b>Bend Points</b></em>' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__BEND_POINTS = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 4;

    /**
     * The feature id for the '<em><b>Parent</b></em>' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__PARENT = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 5;

    /**
     * The feature id for the '<em><b>Outgoing Shape</b></em>' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__OUTGOING_SHAPE = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 6;

    /**
     * The feature id for the '<em><b>Incoming Shape</b></em>' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__INCOMING_SHAPE = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 7;

    /**
     * The feature id for the '<em><b>Outgoing Sections</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__OUTGOING_SECTIONS = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 8;

    /**
     * The feature id for the '<em><b>Incoming Sections</b></em>' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__INCOMING_SECTIONS = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 9;

    /**
     * The feature id for the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION__IDENTIFIER = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 10;

    /**
     * The number of structural features of the '<em>Elk Edge Section</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_EDGE_SECTION_FEATURE_COUNT = EMAP_PROPERTY_HOLDER_FEATURE_COUNT + 11;

    /**
     * The meta object id for the '{@link org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl <em>Elk Property To Value Map Entry</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkPropertyToValueMapEntry()
     * @generated
     */
    int ELK_PROPERTY_TO_VALUE_MAP_ENTRY = 11;

    /**
     * The feature id for the '<em><b>Key</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY = 0;

    /**
     * The feature id for the '<em><b>Value</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE = 1;

    /**
     * The number of structural features of the '<em>Elk Property To Value Map Entry</em>' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     * @ordered
     */
    int ELK_PROPERTY_TO_VALUE_MAP_ENTRY_FEATURE_COUNT = 2;

    /**
     * The meta object id for the '<em>IProperty</em>' data type.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.graph.properties.IProperty
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getIProperty()
     * @generated
     */
    int IPROPERTY = 12;


    /**
     * The meta object id for the '<em>Property Value</em>' data type.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see java.lang.Object
     * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getPropertyValue()
     * @generated
     */
    int PROPERTY_VALUE = 13;


    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.properties.IPropertyHolder <em>IProperty Holder</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>IProperty Holder</em>'.
     * @see org.eclipse.elk.graph.properties.IPropertyHolder
     * @model instanceClass="org.eclipse.elk.graph.properties.IPropertyHolder"
     * @generated
     */
    EClass getIPropertyHolder();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.EMapPropertyHolder <em>EMap Property Holder</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>EMap Property Holder</em>'.
     * @see org.eclipse.elk.graph.EMapPropertyHolder
     * @generated
     */
    EClass getEMapPropertyHolder();

    /**
     * Returns the meta object for the map '{@link org.eclipse.elk.graph.EMapPropertyHolder#getProperties <em>Properties</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the map '<em>Properties</em>'.
     * @see org.eclipse.elk.graph.EMapPropertyHolder#getProperties()
     * @see #getEMapPropertyHolder()
     * @generated
     */
    EReference getEMapPropertyHolder_Properties();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkGraphElement <em>Element</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Element</em>'.
     * @see org.eclipse.elk.graph.ElkGraphElement
     * @generated
     */
    EClass getElkGraphElement();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkGraphElement#getLabels <em>Labels</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Labels</em>'.
     * @see org.eclipse.elk.graph.ElkGraphElement#getLabels()
     * @see #getElkGraphElement()
     * @generated
     */
    EReference getElkGraphElement_Labels();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkGraphElement#getIdentifier <em>Identifier</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Identifier</em>'.
     * @see org.eclipse.elk.graph.ElkGraphElement#getIdentifier()
     * @see #getElkGraphElement()
     * @generated
     */
    EAttribute getElkGraphElement_Identifier();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkShape <em>Elk Shape</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Shape</em>'.
     * @see org.eclipse.elk.graph.ElkShape
     * @generated
     */
    EClass getElkShape();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkShape#getHeight <em>Height</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Height</em>'.
     * @see org.eclipse.elk.graph.ElkShape#getHeight()
     * @see #getElkShape()
     * @generated
     */
    EAttribute getElkShape_Height();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkShape#getWidth <em>Width</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Width</em>'.
     * @see org.eclipse.elk.graph.ElkShape#getWidth()
     * @see #getElkShape()
     * @generated
     */
    EAttribute getElkShape_Width();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkShape#getX <em>X</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>X</em>'.
     * @see org.eclipse.elk.graph.ElkShape#getX()
     * @see #getElkShape()
     * @generated
     */
    EAttribute getElkShape_X();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkShape#getY <em>Y</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Y</em>'.
     * @see org.eclipse.elk.graph.ElkShape#getY()
     * @see #getElkShape()
     * @generated
     */
    EAttribute getElkShape_Y();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkLabel <em>Elk Label</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Label</em>'.
     * @see org.eclipse.elk.graph.ElkLabel
     * @generated
     */
    EClass getElkLabel();

    /**
     * Returns the meta object for the container reference '{@link org.eclipse.elk.graph.ElkLabel#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the container reference '<em>Parent</em>'.
     * @see org.eclipse.elk.graph.ElkLabel#getParent()
     * @see #getElkLabel()
     * @generated
     */
    EReference getElkLabel_Parent();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkLabel#getText <em>Text</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Text</em>'.
     * @see org.eclipse.elk.graph.ElkLabel#getText()
     * @see #getElkLabel()
     * @generated
     */
    EAttribute getElkLabel_Text();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkConnectableShape <em>Elk Connectable Shape</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Connectable Shape</em>'.
     * @see org.eclipse.elk.graph.ElkConnectableShape
     * @generated
     */
    EClass getElkConnectableShape();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkConnectableShape#getOutgoingEdges <em>Outgoing Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Outgoing Edges</em>'.
     * @see org.eclipse.elk.graph.ElkConnectableShape#getOutgoingEdges()
     * @see #getElkConnectableShape()
     * @generated
     */
    EReference getElkConnectableShape_OutgoingEdges();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkConnectableShape#getIncomingEdges <em>Incoming Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Incoming Edges</em>'.
     * @see org.eclipse.elk.graph.ElkConnectableShape#getIncomingEdges()
     * @see #getElkConnectableShape()
     * @generated
     */
    EReference getElkConnectableShape_IncomingEdges();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkNode <em>Elk Node</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Node</em>'.
     * @see org.eclipse.elk.graph.ElkNode
     * @generated
     */
    EClass getElkNode();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkNode#getPorts <em>Ports</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Ports</em>'.
     * @see org.eclipse.elk.graph.ElkNode#getPorts()
     * @see #getElkNode()
     * @generated
     */
    EReference getElkNode_Ports();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkNode#getChildren <em>Children</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Children</em>'.
     * @see org.eclipse.elk.graph.ElkNode#getChildren()
     * @see #getElkNode()
     * @generated
     */
    EReference getElkNode_Children();

    /**
     * Returns the meta object for the container reference '{@link org.eclipse.elk.graph.ElkNode#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the container reference '<em>Parent</em>'.
     * @see org.eclipse.elk.graph.ElkNode#getParent()
     * @see #getElkNode()
     * @generated
     */
    EReference getElkNode_Parent();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkNode#getContainedEdges <em>Contained Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Contained Edges</em>'.
     * @see org.eclipse.elk.graph.ElkNode#getContainedEdges()
     * @see #getElkNode()
     * @generated
     */
    EReference getElkNode_ContainedEdges();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkNode#isHierarchical <em>Hierarchical</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Hierarchical</em>'.
     * @see org.eclipse.elk.graph.ElkNode#isHierarchical()
     * @see #getElkNode()
     * @generated
     */
    EAttribute getElkNode_Hierarchical();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkPort <em>Elk Port</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Port</em>'.
     * @see org.eclipse.elk.graph.ElkPort
     * @generated
     */
    EClass getElkPort();

    /**
     * Returns the meta object for the container reference '{@link org.eclipse.elk.graph.ElkPort#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the container reference '<em>Parent</em>'.
     * @see org.eclipse.elk.graph.ElkPort#getParent()
     * @see #getElkPort()
     * @generated
     */
    EReference getElkPort_Parent();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkEdge <em>Elk Edge</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Edge</em>'.
     * @see org.eclipse.elk.graph.ElkEdge
     * @generated
     */
    EClass getElkEdge();

    /**
     * Returns the meta object for the container reference '{@link org.eclipse.elk.graph.ElkEdge#getContainingNode <em>Containing Node</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the container reference '<em>Containing Node</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#getContainingNode()
     * @see #getElkEdge()
     * @generated
     */
    EReference getElkEdge_ContainingNode();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkEdge#getSources <em>Sources</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Sources</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#getSources()
     * @see #getElkEdge()
     * @generated
     */
    EReference getElkEdge_Sources();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkEdge#getTargets <em>Targets</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Targets</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#getTargets()
     * @see #getElkEdge()
     * @generated
     */
    EReference getElkEdge_Targets();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkEdge#getSections <em>Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Sections</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#getSections()
     * @see #getElkEdge()
     * @generated
     */
    EReference getElkEdge_Sections();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdge#isHyperedge <em>Hyperedge</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Hyperedge</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#isHyperedge()
     * @see #getElkEdge()
     * @generated
     */
    EAttribute getElkEdge_Hyperedge();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdge#isHierarchical <em>Hierarchical</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Hierarchical</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#isHierarchical()
     * @see #getElkEdge()
     * @generated
     */
    EAttribute getElkEdge_Hierarchical();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdge#isSelfloop <em>Selfloop</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Selfloop</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#isSelfloop()
     * @see #getElkEdge()
     * @generated
     */
    EAttribute getElkEdge_Selfloop();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdge#isConnected <em>Connected</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Connected</em>'.
     * @see org.eclipse.elk.graph.ElkEdge#isConnected()
     * @see #getElkEdge()
     * @generated
     */
    EAttribute getElkEdge_Connected();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkBendPoint <em>Elk Bend Point</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Bend Point</em>'.
     * @see org.eclipse.elk.graph.ElkBendPoint
     * @generated
     */
    EClass getElkBendPoint();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkBendPoint#getX <em>X</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>X</em>'.
     * @see org.eclipse.elk.graph.ElkBendPoint#getX()
     * @see #getElkBendPoint()
     * @generated
     */
    EAttribute getElkBendPoint_X();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkBendPoint#getY <em>Y</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Y</em>'.
     * @see org.eclipse.elk.graph.ElkBendPoint#getY()
     * @see #getElkBendPoint()
     * @generated
     */
    EAttribute getElkBendPoint_Y();

    /**
     * Returns the meta object for class '{@link org.eclipse.elk.graph.ElkEdgeSection <em>Elk Edge Section</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Edge Section</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection
     * @generated
     */
    EClass getElkEdgeSection();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdgeSection#getStartX <em>Start X</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Start X</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getStartX()
     * @see #getElkEdgeSection()
     * @generated
     */
    EAttribute getElkEdgeSection_StartX();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdgeSection#getStartY <em>Start Y</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Start Y</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getStartY()
     * @see #getElkEdgeSection()
     * @generated
     */
    EAttribute getElkEdgeSection_StartY();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdgeSection#getEndX <em>End X</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>End X</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getEndX()
     * @see #getElkEdgeSection()
     * @generated
     */
    EAttribute getElkEdgeSection_EndX();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdgeSection#getEndY <em>End Y</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>End Y</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getEndY()
     * @see #getElkEdgeSection()
     * @generated
     */
    EAttribute getElkEdgeSection_EndY();

    /**
     * Returns the meta object for the containment reference list '{@link org.eclipse.elk.graph.ElkEdgeSection#getBendPoints <em>Bend Points</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the containment reference list '<em>Bend Points</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getBendPoints()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_BendPoints();

    /**
     * Returns the meta object for the container reference '{@link org.eclipse.elk.graph.ElkEdgeSection#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the container reference '<em>Parent</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getParent()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_Parent();

    /**
     * Returns the meta object for the reference '{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingShape <em>Outgoing Shape</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference '<em>Outgoing Shape</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getOutgoingShape()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_OutgoingShape();

    /**
     * Returns the meta object for the reference '{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingShape <em>Incoming Shape</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference '<em>Incoming Shape</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getIncomingShape()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_IncomingShape();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingSections <em>Outgoing Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Outgoing Sections</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getOutgoingSections()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_OutgoingSections();

    /**
     * Returns the meta object for the reference list '{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingSections <em>Incoming Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the reference list '<em>Incoming Sections</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getIncomingSections()
     * @see #getElkEdgeSection()
     * @generated
     */
    EReference getElkEdgeSection_IncomingSections();

    /**
     * Returns the meta object for the attribute '{@link org.eclipse.elk.graph.ElkEdgeSection#getIdentifier <em>Identifier</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Identifier</em>'.
     * @see org.eclipse.elk.graph.ElkEdgeSection#getIdentifier()
     * @see #getElkEdgeSection()
     * @generated
     */
    EAttribute getElkEdgeSection_Identifier();

    /**
     * Returns the meta object for class '{@link java.util.Map.Entry <em>Elk Property To Value Map Entry</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for class '<em>Elk Property To Value Map Entry</em>'.
     * @see java.util.Map.Entry
     * @model keyDataType="org.eclipse.elk.graph.IProperty&lt;?&gt;"
     *        valueDataType="org.eclipse.elk.graph.PropertyValue"
     * @generated
     */
    EClass getElkPropertyToValueMapEntry();

    /**
     * Returns the meta object for the attribute '{@link java.util.Map.Entry <em>Key</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Key</em>'.
     * @see java.util.Map.Entry
     * @see #getElkPropertyToValueMapEntry()
     * @generated
     */
    EAttribute getElkPropertyToValueMapEntry_Key();

    /**
     * Returns the meta object for the attribute '{@link java.util.Map.Entry <em>Value</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for the attribute '<em>Value</em>'.
     * @see java.util.Map.Entry
     * @see #getElkPropertyToValueMapEntry()
     * @generated
     */
    EAttribute getElkPropertyToValueMapEntry_Value();

    /**
     * Returns the meta object for data type '{@link org.eclipse.elk.graph.properties.IProperty <em>IProperty</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for data type '<em>IProperty</em>'.
     * @see org.eclipse.elk.graph.properties.IProperty
     * @model instanceClass="org.eclipse.elk.graph.properties.IProperty" typeParameters="T"
     * @generated
     */
    EDataType getIProperty();

    /**
     * Returns the meta object for data type '{@link java.lang.Object <em>Property Value</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the meta object for data type '<em>Property Value</em>'.
     * @see java.lang.Object
     * @model instanceClass="java.lang.Object"
     * @generated
     */
    EDataType getPropertyValue();

    /**
     * Returns the factory that creates the instances of the model.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the factory that creates the instances of the model.
     * @generated
     */
    ElkGraphFactory getElkGraphFactory();

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
    interface Literals {
        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.properties.IPropertyHolder <em>IProperty Holder</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.properties.IPropertyHolder
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getIPropertyHolder()
         * @generated
         */
        EClass IPROPERTY_HOLDER = eINSTANCE.getIPropertyHolder();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.EMapPropertyHolderImpl <em>EMap Property Holder</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.EMapPropertyHolderImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getEMapPropertyHolder()
         * @generated
         */
        EClass EMAP_PROPERTY_HOLDER = eINSTANCE.getEMapPropertyHolder();

        /**
         * The meta object literal for the '<em><b>Properties</b></em>' map feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference EMAP_PROPERTY_HOLDER__PROPERTIES = eINSTANCE.getEMapPropertyHolder_Properties();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkGraphElementImpl <em>Element</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkGraphElementImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkGraphElement()
         * @generated
         */
        EClass ELK_GRAPH_ELEMENT = eINSTANCE.getElkGraphElement();

        /**
         * The meta object literal for the '<em><b>Labels</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_GRAPH_ELEMENT__LABELS = eINSTANCE.getElkGraphElement_Labels();

        /**
         * The meta object literal for the '<em><b>Identifier</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_GRAPH_ELEMENT__IDENTIFIER = eINSTANCE.getElkGraphElement_Identifier();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkShapeImpl <em>Elk Shape</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkShapeImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkShape()
         * @generated
         */
        EClass ELK_SHAPE = eINSTANCE.getElkShape();

        /**
         * The meta object literal for the '<em><b>Height</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_SHAPE__HEIGHT = eINSTANCE.getElkShape_Height();

        /**
         * The meta object literal for the '<em><b>Width</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_SHAPE__WIDTH = eINSTANCE.getElkShape_Width();

        /**
         * The meta object literal for the '<em><b>X</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_SHAPE__X = eINSTANCE.getElkShape_X();

        /**
         * The meta object literal for the '<em><b>Y</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_SHAPE__Y = eINSTANCE.getElkShape_Y();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkLabelImpl <em>Elk Label</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkLabelImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkLabel()
         * @generated
         */
        EClass ELK_LABEL = eINSTANCE.getElkLabel();

        /**
         * The meta object literal for the '<em><b>Parent</b></em>' container reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_LABEL__PARENT = eINSTANCE.getElkLabel_Parent();

        /**
         * The meta object literal for the '<em><b>Text</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_LABEL__TEXT = eINSTANCE.getElkLabel_Text();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkConnectableShapeImpl <em>Elk Connectable Shape</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkConnectableShapeImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkConnectableShape()
         * @generated
         */
        EClass ELK_CONNECTABLE_SHAPE = eINSTANCE.getElkConnectableShape();

        /**
         * The meta object literal for the '<em><b>Outgoing Edges</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES = eINSTANCE.getElkConnectableShape_OutgoingEdges();

        /**
         * The meta object literal for the '<em><b>Incoming Edges</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_CONNECTABLE_SHAPE__INCOMING_EDGES = eINSTANCE.getElkConnectableShape_IncomingEdges();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkNodeImpl <em>Elk Node</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkNodeImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkNode()
         * @generated
         */
        EClass ELK_NODE = eINSTANCE.getElkNode();

        /**
         * The meta object literal for the '<em><b>Ports</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_NODE__PORTS = eINSTANCE.getElkNode_Ports();

        /**
         * The meta object literal for the '<em><b>Children</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_NODE__CHILDREN = eINSTANCE.getElkNode_Children();

        /**
         * The meta object literal for the '<em><b>Parent</b></em>' container reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_NODE__PARENT = eINSTANCE.getElkNode_Parent();

        /**
         * The meta object literal for the '<em><b>Contained Edges</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_NODE__CONTAINED_EDGES = eINSTANCE.getElkNode_ContainedEdges();

        /**
         * The meta object literal for the '<em><b>Hierarchical</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_NODE__HIERARCHICAL = eINSTANCE.getElkNode_Hierarchical();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkPortImpl <em>Elk Port</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkPortImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkPort()
         * @generated
         */
        EClass ELK_PORT = eINSTANCE.getElkPort();

        /**
         * The meta object literal for the '<em><b>Parent</b></em>' container reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_PORT__PARENT = eINSTANCE.getElkPort_Parent();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkEdgeImpl <em>Elk Edge</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkEdgeImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkEdge()
         * @generated
         */
        EClass ELK_EDGE = eINSTANCE.getElkEdge();

        /**
         * The meta object literal for the '<em><b>Containing Node</b></em>' container reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE__CONTAINING_NODE = eINSTANCE.getElkEdge_ContainingNode();

        /**
         * The meta object literal for the '<em><b>Sources</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE__SOURCES = eINSTANCE.getElkEdge_Sources();

        /**
         * The meta object literal for the '<em><b>Targets</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE__TARGETS = eINSTANCE.getElkEdge_Targets();

        /**
         * The meta object literal for the '<em><b>Sections</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE__SECTIONS = eINSTANCE.getElkEdge_Sections();

        /**
         * The meta object literal for the '<em><b>Hyperedge</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE__HYPEREDGE = eINSTANCE.getElkEdge_Hyperedge();

        /**
         * The meta object literal for the '<em><b>Hierarchical</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE__HIERARCHICAL = eINSTANCE.getElkEdge_Hierarchical();

        /**
         * The meta object literal for the '<em><b>Selfloop</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE__SELFLOOP = eINSTANCE.getElkEdge_Selfloop();

        /**
         * The meta object literal for the '<em><b>Connected</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE__CONNECTED = eINSTANCE.getElkEdge_Connected();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkBendPointImpl <em>Elk Bend Point</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkBendPointImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkBendPoint()
         * @generated
         */
        EClass ELK_BEND_POINT = eINSTANCE.getElkBendPoint();

        /**
         * The meta object literal for the '<em><b>X</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_BEND_POINT__X = eINSTANCE.getElkBendPoint_X();

        /**
         * The meta object literal for the '<em><b>Y</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_BEND_POINT__Y = eINSTANCE.getElkBendPoint_Y();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl <em>Elk Edge Section</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkEdgeSectionImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkEdgeSection()
         * @generated
         */
        EClass ELK_EDGE_SECTION = eINSTANCE.getElkEdgeSection();

        /**
         * The meta object literal for the '<em><b>Start X</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE_SECTION__START_X = eINSTANCE.getElkEdgeSection_StartX();

        /**
         * The meta object literal for the '<em><b>Start Y</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE_SECTION__START_Y = eINSTANCE.getElkEdgeSection_StartY();

        /**
         * The meta object literal for the '<em><b>End X</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE_SECTION__END_X = eINSTANCE.getElkEdgeSection_EndX();

        /**
         * The meta object literal for the '<em><b>End Y</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE_SECTION__END_Y = eINSTANCE.getElkEdgeSection_EndY();

        /**
         * The meta object literal for the '<em><b>Bend Points</b></em>' containment reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__BEND_POINTS = eINSTANCE.getElkEdgeSection_BendPoints();

        /**
         * The meta object literal for the '<em><b>Parent</b></em>' container reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__PARENT = eINSTANCE.getElkEdgeSection_Parent();

        /**
         * The meta object literal for the '<em><b>Outgoing Shape</b></em>' reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__OUTGOING_SHAPE = eINSTANCE.getElkEdgeSection_OutgoingShape();

        /**
         * The meta object literal for the '<em><b>Incoming Shape</b></em>' reference feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__INCOMING_SHAPE = eINSTANCE.getElkEdgeSection_IncomingShape();

        /**
         * The meta object literal for the '<em><b>Outgoing Sections</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__OUTGOING_SECTIONS = eINSTANCE.getElkEdgeSection_OutgoingSections();

        /**
         * The meta object literal for the '<em><b>Incoming Sections</b></em>' reference list feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EReference ELK_EDGE_SECTION__INCOMING_SECTIONS = eINSTANCE.getElkEdgeSection_IncomingSections();

        /**
         * The meta object literal for the '<em><b>Identifier</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_EDGE_SECTION__IDENTIFIER = eINSTANCE.getElkEdgeSection_Identifier();

        /**
         * The meta object literal for the '{@link org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl <em>Elk Property To Value Map Entry</em>}' class.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getElkPropertyToValueMapEntry()
         * @generated
         */
        EClass ELK_PROPERTY_TO_VALUE_MAP_ENTRY = eINSTANCE.getElkPropertyToValueMapEntry();

        /**
         * The meta object literal for the '<em><b>Key</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY = eINSTANCE.getElkPropertyToValueMapEntry_Key();

        /**
         * The meta object literal for the '<em><b>Value</b></em>' attribute feature.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @generated
         */
        EAttribute ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE = eINSTANCE.getElkPropertyToValueMapEntry_Value();

        /**
         * The meta object literal for the '<em>IProperty</em>' data type.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see org.eclipse.elk.graph.properties.IProperty
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getIProperty()
         * @generated
         */
        EDataType IPROPERTY = eINSTANCE.getIProperty();

        /**
         * The meta object literal for the '<em>Property Value</em>' data type.
         * <!-- begin-user-doc -->
         * <!-- end-user-doc -->
         * @see java.lang.Object
         * @see org.eclipse.elk.graph.impl.ElkGraphPackageImpl#getPropertyValue()
         * @generated
         */
        EDataType PROPERTY_VALUE = eINSTANCE.getPropertyValue();

    }

} //ElkGraphPackage
