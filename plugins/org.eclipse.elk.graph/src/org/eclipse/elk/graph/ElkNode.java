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

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Node</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A node in the graph. Edges can be connected to the node directly or through one of its ports.
 * 
 * <p>All nodes except one must have an assigned parent node. The node that does not have a parent node is the graph's root node and represents the graph itself. There can only be one root node for each graph. The parent-child relationship of nodes induces hierarchy in nested graphs: a node's children constitute the graph contained in and represented by that node.</p>
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkNode#getPorts <em>Ports</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkNode#getChildren <em>Children</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkNode#getParent <em>Parent</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkNode#getContainedEdges <em>Contained Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkNode#isHierarchical <em>Hierarchical</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode()
 * @model
 * @generated
 */
public interface ElkNode extends ElkConnectableShape {
    /**
     * Returns the value of the '<em><b>Ports</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkPort}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkPort#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The node's list of ports.
     * 
     * <p>Adding or removing a port to/from this list automatically sets its parent node.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Ports</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode_Ports()
     * @see org.eclipse.elk.graph.ElkPort#getParent
     * @model opposite="parent" containment="true"
     * @generated
     */
    EList<ElkPort> getPorts();

    /**
     * Returns the value of the '<em><b>Children</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkNode}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkNode#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Child nodes contained in this node. If the node contains at least one child node, the node is a <em>hierarchical node</em>.
     * 
     * <p>Adding or removing a node to/from this list automatically sets its parent node.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Children</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode_Children()
     * @see org.eclipse.elk.graph.ElkNode#getParent
     * @model opposite="parent" containment="true"
     * @generated
     */
    EList<ElkNode> getChildren();

    /**
     * Returns the value of the '<em><b>Parent</b></em>' container reference.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkNode#getChildren <em>Children</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The node's parent node, if any. If the node has a parent, its location is interpreted relative to the parent node's origin. If it does not, the meaning of its location is not defined. A node without a parent node is called the <em>root node</em>.
     * 
     * <p>Setting the node's parent node automatically updates the parent node's list of child nodes.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Parent</em>' container reference.
     * @see #setParent(ElkNode)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode_Parent()
     * @see org.eclipse.elk.graph.ElkNode#getChildren
     * @model opposite="children" transient="false"
     * @generated
     */
    ElkNode getParent();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkNode#getParent <em>Parent</em>}' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Parent</em>' container reference.
     * @see #getParent()
     * @generated
     */
    void setParent(ElkNode value);

    /**
     * Returns the value of the '<em><b>Contained Edges</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdge}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdge#getContainingNode <em>Containing Node</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The edges contained in this node. Note that containment does not imply connection: an edge contained in this list is not necessarily connected to this node. The coordinates of the edge's source points, target points, and end points are interpreted relative to this node's origin. See the documentation for {@link ElkEdge} for more details on which node an edge should be contained in.
     * 
     * <p>Adding or removing an edge to/from this list automatically sets its containing node.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Contained Edges</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode_ContainedEdges()
     * @see org.eclipse.elk.graph.ElkEdge#getContainingNode
     * @model opposite="containingNode" containment="true"
     * @generated
     */
    EList<ElkEdge> getContainedEdges();

    /**
     * Returns the value of the '<em><b>Hierarchical</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Whether or not this node is considered to be a hierarchical node. It is a hierarchical node if its list of child nodes is not empty.
     * 
     * <p>The value of this attribute is computed dynamically and not persistent. So, no worries if you can't find it anywhere.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Hierarchical</em>' attribute.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkNode_Hierarchical()
     * @model transient="true" changeable="false" volatile="true" derived="true"
     * @generated
     */
    boolean isHierarchical();

} // ElkNode
