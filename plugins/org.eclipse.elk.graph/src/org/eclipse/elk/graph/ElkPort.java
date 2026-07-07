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


/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Port</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A port represents an explicit point through which to connect to a node. Different ports of a node will usually have different associated meanings, much like different method parameters. Each port belongs to the node it is contained in.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkPort#getParent <em>Parent</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkPort()
 * @model
 * @generated
 */
public interface ElkPort extends ElkConnectableShape {
    /**
     * Returns the value of the '<em><b>Parent</b></em>' container reference.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkNode#getPorts <em>Ports</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The node the port belongs to.
     * 
     * <p>Setting the parent node automatically update's the node's list of ports.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Parent</em>' container reference.
     * @see #setParent(ElkNode)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkPort_Parent()
     * @see org.eclipse.elk.graph.ElkNode#getPorts
     * @model opposite="ports" transient="false"
     * @generated
     */
    ElkNode getParent();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkPort#getParent <em>Parent</em>}' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Parent</em>' container reference.
     * @see #getParent()
     * @generated
     */
    void setParent(ElkNode value);

} // ElkPort
