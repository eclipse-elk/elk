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
 * A representation of the model object '<em><b>Elk Edge</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * An edge connects one or more source elements ({@link ElkConnectableShape}s) with one or more target elements ({@link ElkConnectableShape}s). If an edge connects at most one source with at most one target, it is called a <em>regular edge</em> (although it is usually simply called an <em>edge</em>). If an edge has more than a single source or more than a single target, it is called a <em>hyperedge</em>. If all of the edge's sources and targets have the same parent node, it is a <em>simple edge</em>; otherwise, it is called a <em>hierarchical edge</em>.
 * 
 * Each edge must be assigned to a containing node. The containing node defines the point where it is hooked into the graph's object hierarchy, which is important for serializing the graph. The containing node's origin is the point the edge's source, target, and bend points are relative to. As a rule of thumb, edges should always be contained in the lowest common representing node of the graphs of all elements it connects, with one exception: if an edge connects a node with one of its descendants, that node should be the edge's containing node. The {@link GraphUtil} class contains methods that help with finding the correct containing node.
 * 
 * The routing of an edge is specified by the {@link ElkEdgeSection} objects it contains. If the edge is a regular edge (as opposed to a hyperedge), it contains at most a single {@link ElkEdgeSection} which specifies a single source point, a single end point, and an arbitrary number of bend points. If the edge is a hyperedge, it contains at least one {@link ElkEdgeSection} for each of its sources and targets (<em>outer edge sections</em>) as well as an arbitrary number of {@link ElkEdgeSection} objects to connect the outer sections (<em>inner edge sections</em>).
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#getContainingNode <em>Containing Node</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#getSources <em>Sources</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#getTargets <em>Targets</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#getSections <em>Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#isHyperedge <em>Hyperedge</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#isHierarchical <em>Hierarchical</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#isSelfloop <em>Selfloop</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdge#isConnected <em>Connected</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge()
 * @model
 * @generated
 */
public interface ElkEdge extends ElkGraphElement {
    /**
     * Returns the value of the '<em><b>Containing Node</b></em>' container reference.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkNode#getContainedEdges <em>Contained Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The node the edge is contained in. This has nothing to do with what the edge connects to.
     * 
     * <p>Setting the contained node automatically update's the node's list of contained edges.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Containing Node</em>' container reference.
     * @see #setContainingNode(ElkNode)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_ContainingNode()
     * @see org.eclipse.elk.graph.ElkNode#getContainedEdges
     * @model opposite="containedEdges" transient="false"
     * @generated
     */
    ElkNode getContainingNode();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdge#getContainingNode <em>Containing Node</em>}' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Containing Node</em>' container reference.
     * @see #getContainingNode()
     * @generated
     */
    void setContainingNode(ElkNode value);

    /**
     * Returns the value of the '<em><b>Sources</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkConnectableShape}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkConnectableShape#getOutgoingEdges <em>Outgoing Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The edge's list of source elements. For an edge to make sense, it should have at least one source element.
     * 
     * <p>Adding or removing an element to/from this list automatically updates its list of outgoing edges.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Sources</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Sources()
     * @see org.eclipse.elk.graph.ElkConnectableShape#getOutgoingEdges
     * @model opposite="outgoingEdges"
     * @generated
     */
    EList<ElkConnectableShape> getSources();

    /**
     * Returns the value of the '<em><b>Targets</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkConnectableShape}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkConnectableShape#getIncomingEdges <em>Incoming Edges</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The edge's list of target elements. For an edge to make sense, it should have at least one target element.
     * 
     * <p>Adding or removing an element to/from this list automatically updates its list of incoming edges.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Targets</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Targets()
     * @see org.eclipse.elk.graph.ElkConnectableShape#getIncomingEdges
     * @model opposite="incomingEdges"
     * @generated
     */
    EList<ElkConnectableShape> getTargets();

    /**
     * Returns the value of the '<em><b>Sections</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdgeSection}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdgeSection#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * All edge sections that define the routing of this edge.
     * 
     * <p>Adding or removing an edge section to/from this list automatically updates its parent edge.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Sections</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Sections()
     * @see org.eclipse.elk.graph.ElkEdgeSection#getParent
     * @model opposite="parent" containment="true"
     * @generated
     */
    EList<ElkEdgeSection> getSections();

    /**
     * Returns the value of the '<em><b>Hyperedge</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Whether this edge is a hyperedge or not. It is considered to not be a hyperedge if it has at most one source and at most one target. The result is only meaningful if the node has at least one source and one target.
     * 
     * <p>The value of this attribute is computed dynamically and not persistent. So, no worries if you can't find it anywhere.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Hyperedge</em>' attribute.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Hyperedge()
     * @model transient="true" changeable="false" volatile="true" derived="true"
     * @generated
     */
    boolean isHyperedge();

    /**
     * Returns the value of the '<em><b>Hierarchical</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Whether the edge is a hierarchical edge or not. It is hierarchical if its sources and targets do not all have the same parent node.
     * 
     * <p>The value of this attribute is computed dynamically and not persistent. So, no worries if you can't find it anywhere.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Hierarchical</em>' attribute.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Hierarchical()
     * @model transient="true" changeable="false" volatile="true" derived="true"
     * @generated
     */
    boolean isHierarchical();

    /**
     * Returns the value of the '<em><b>Selfloop</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Whether the edge is a self loop or not. It is a self loop if all of its sources and targets are the same node (or its ports). This is only meaningful if the edge connects at least two end points.
     * 
     * <p>The value of this attribute is computed dynamically and not persistent. So, no worries if you can't find it anywhere.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Selfloop</em>' attribute.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Selfloop()
     * @model transient="true" changeable="false" volatile="true" derived="true"
     * @generated
     */
    boolean isSelfloop();

    /**
     * Returns the value of the '<em><b>Connected</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Whether the edge has at least one source and one target.
     * 
     * <p>The value of this attribute is computed dynamically and not persistent. So, no worries if you can't find it anywhere.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Connected</em>' attribute.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdge_Connected()
     * @model transient="true" changeable="false" volatile="true" derived="true"
     * @generated
     */
    boolean isConnected();

} // ElkEdge
