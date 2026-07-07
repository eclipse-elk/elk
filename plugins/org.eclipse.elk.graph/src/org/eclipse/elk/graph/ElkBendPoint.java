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

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Bend Point</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A bend point of an {@link ElkEdgeSection}. The coordinates of a bend point are always relative to the origin of the containing node of the edge the bend point belongs to.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkBendPoint#getX <em>X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkBendPoint#getY <em>Y</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkBendPoint()
 * @model
 * @generated
 */
public interface ElkBendPoint extends EObject {
    /**
     * Returns the value of the '<em><b>X</b></em>' attribute.
     * The default value is <code>"0.0"</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The bend point's x coordinate, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>X</em>' attribute.
     * @see #setX(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkBendPoint_X()
     * @model default="0.0" required="true"
     * @generated
     */
    double getX();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkBendPoint#getX <em>X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>X</em>' attribute.
     * @see #getX()
     * @generated
     */
    void setX(double value);

    /**
     * Returns the value of the '<em><b>Y</b></em>' attribute.
     * The default value is <code>"0.0"</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The bend point's y coordinate, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Y</em>' attribute.
     * @see #setY(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkBendPoint_Y()
     * @model default="0.0" required="true"
     * @generated
     */
    double getY();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkBendPoint#getY <em>Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Y</em>' attribute.
     * @see #getY()
     * @generated
     */
    void setY(double value);

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Sets the bend point's x and y coordinates simultaneously by calling their respective set methods.
     * <!-- end-model-doc -->
     * @model annotation="http://www.eclipse.org/emf/2002/GenModel body='setX(x);\nsetY(y);'"
     * @generated
     */
    void set(double x, double y);

} // ElkBendPoint
