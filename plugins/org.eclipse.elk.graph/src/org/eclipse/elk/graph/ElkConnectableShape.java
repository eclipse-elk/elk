/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph;

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Connectable Shape</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A graph element that can be the end point of an edge.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkConnectableShape#getOutgoingEdges <em>Outgoing Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkConnectableShape#getIncomingEdges <em>Incoming Edges</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkConnectableShape()
 * @model abstract="true"
 * @generated
 */
public interface ElkConnectableShape extends ElkShape {
    /**
     * Returns the value of the '<em><b>Outgoing Edges</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdge}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdge#getSources <em>Sources</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * List of edges that leave this connectable shape.
     * 
     * <p>Adding or removing an edge to/from this list automatically updates its list of sources.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Outgoing Edges</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkConnectableShape_OutgoingEdges()
     * @see org.eclipse.elk.graph.ElkEdge#getSources
     * @model opposite="sources"
     * @generated
     */
    EList<ElkEdge> getOutgoingEdges();

    /**
     * Returns the value of the '<em><b>Incoming Edges</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdge}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdge#getTargets <em>Targets</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * List of edges that go into this connectable shape.
     * 
     * <p>Adding or removing an edge to/from this list automatically updates its list of targets.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Incoming Edges</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkConnectableShape_IncomingEdges()
     * @see org.eclipse.elk.graph.ElkEdge#getTargets
     * @model opposite="targets"
     * @generated
     */
    EList<ElkEdge> getIncomingEdges();

} // ElkConnectableShape
