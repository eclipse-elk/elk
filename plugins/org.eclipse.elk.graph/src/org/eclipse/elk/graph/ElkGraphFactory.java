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

import org.eclipse.emf.ecore.EFactory;

/**
 * <!-- begin-user-doc -->
 * The <b>Factory</b> for the model.
 * It provides a create method for each non-abstract class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.graph.ElkGraphPackage
 * @generated
 */
public interface ElkGraphFactory extends EFactory {
    /**
     * The singleton instance of the factory.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    ElkGraphFactory eINSTANCE = org.eclipse.elk.graph.impl.ElkGraphFactoryImpl.init();

    /**
     * Returns a new object of class '<em>Elk Label</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Label</em>'.
     * @generated
     */
    ElkLabel createElkLabel();

    /**
     * Returns a new object of class '<em>Elk Node</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Node</em>'.
     * @generated
     */
    ElkNode createElkNode();

    /**
     * Returns a new object of class '<em>Elk Port</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Port</em>'.
     * @generated
     */
    ElkPort createElkPort();

    /**
     * Returns a new object of class '<em>Elk Edge</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Edge</em>'.
     * @generated
     */
    ElkEdge createElkEdge();

    /**
     * Returns a new object of class '<em>Elk Bend Point</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Bend Point</em>'.
     * @generated
     */
    ElkBendPoint createElkBendPoint();

    /**
     * Returns a new object of class '<em>Elk Edge Section</em>'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return a new object of class '<em>Elk Edge Section</em>'.
     * @generated
     */
    ElkEdgeSection createElkEdgeSection();

    /**
     * Returns the package supported by this factory.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the package supported by this factory.
     * @generated
     */
    ElkGraphPackage getElkGraphPackage();

} //ElkGraphFactory
