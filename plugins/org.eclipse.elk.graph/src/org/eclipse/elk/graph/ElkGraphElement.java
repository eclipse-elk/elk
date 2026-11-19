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
 * A representation of the model object '<em><b>Element</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * This is the superclass of all elements of a graph such as nodes, edges, ports, and labels. Each element can have an arbitrary number of labels attached to it. A graph element can also hold properties that, for instance, influence how it is treated by layout algorithms. Finally, each graph element can have an arbitrary number of {@link ElkGraphData} objects associated with it to further annotate it with more specific information.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkGraphElement#getLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkGraphElement#getIdentifier <em>Identifier</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkGraphElement()
 * @model abstract="true"
 * @generated
 */
public interface ElkGraphElement extends EMapPropertyHolder {
    /**
     * Returns the value of the '<em><b>Labels</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkLabel}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkLabel#getParent <em>Parent</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Labels associated with this graph element.
     * 
     * <p>Adding or removing a label to/from this list automatically updates its parent element.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Labels</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkGraphElement_Labels()
     * @see org.eclipse.elk.graph.ElkLabel#getParent
     * @model opposite="parent" containment="true"
     * @generated
     */
    EList<ElkLabel> getLabels();

    /**
     * Returns the value of the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * An optional String identifier for this graph element. Can be used as an ID for defining graphs in Xtext-based languages.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Identifier</em>' attribute.
     * @see #setIdentifier(String)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkGraphElement_Identifier()
     * @model
     * @generated
     */
    String getIdentifier();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkGraphElement#getIdentifier <em>Identifier</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Identifier</em>' attribute.
     * @see #getIdentifier()
     * @generated
     */
    void setIdentifier(String value);

} // ElkGraphElement
