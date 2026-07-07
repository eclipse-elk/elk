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
 * A representation of the model object '<em><b>Elk Label</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * Labels are used to associate graph elements with text to be displayed in a diagram. The element the label annotates is its parent element.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkLabel#getParent <em>Parent</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkLabel#getText <em>Text</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkLabel()
 * @model
 * @generated
 */
public interface ElkLabel extends ElkShape {
    /**
     * Returns the value of the '<em><b>Parent</b></em>' container reference.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkGraphElement#getLabels <em>Labels</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Graph element the label annotates.
     * 
     * <p>Setting the parent element automatically updates its list of labels.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Parent</em>' container reference.
     * @see #setParent(ElkGraphElement)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkLabel_Parent()
     * @see org.eclipse.elk.graph.ElkGraphElement#getLabels
     * @model opposite="labels" transient="false"
     * @generated
     */
    ElkGraphElement getParent();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkLabel#getParent <em>Parent</em>}' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Parent</em>' container reference.
     * @see #getParent()
     * @generated
     */
    void setParent(ElkGraphElement value);

    /**
     * Returns the value of the '<em><b>Text</b></em>' attribute.
     * The default value is <code>""</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The label's text.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Text</em>' attribute.
     * @see #setText(String)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkLabel_Text()
     * @model default=""
     * @generated
     */
    String getText();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkLabel#getText <em>Text</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Text</em>' attribute.
     * @see #getText()
     * @generated
     */
    void setText(String value);

} // ElkLabel
