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


/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Shape</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A shape is a graph element whose placement and extend can be described by x and y coordinates as well as a width and a height. The coordinates of a shape describe the position of its top left corner, relative to the origin of its parent element. The width and height of a shape describe the extend of its rectangular bounding box.
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkShape#getHeight <em>Height</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkShape#getWidth <em>Width</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkShape#getX <em>X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkShape#getY <em>Y</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkShape()
 * @model abstract="true"
 * @generated
 */
public interface ElkShape extends ElkGraphElement {
    /**
     * Returns the value of the '<em><b>Height</b></em>' attribute.
     * The default value is <code>"0.0"</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Height of the shape's rectangular bounding box.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Height</em>' attribute.
     * @see #setHeight(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkShape_Height()
     * @model default="0.0" required="true"
     * @generated
     */
    double getHeight();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkShape#getHeight <em>Height</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Height</em>' attribute.
     * @see #getHeight()
     * @generated
     */
    void setHeight(double value);

    /**
     * Returns the value of the '<em><b>Width</b></em>' attribute.
     * The default value is <code>"0.0"</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Height of the shape's rectangular bounding box.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Width</em>' attribute.
     * @see #setWidth(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkShape_Width()
     * @model default="0.0" required="true"
     * @generated
     */
    double getWidth();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkShape#getWidth <em>Width</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Width</em>' attribute.
     * @see #getWidth()
     * @generated
     */
    void setWidth(double value);

    /**
     * Returns the value of the '<em><b>X</b></em>' attribute.
     * The default value is <code>"0.0"</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * X coordinate of the shape's top left corner, relative to the origin of its parent object.
     * <!-- end-model-doc -->
     * @return the value of the '<em>X</em>' attribute.
     * @see #setX(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkShape_X()
     * @model default="0.0" required="true"
     * @generated
     */
    double getX();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkShape#getX <em>X</em>}' attribute.
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
     * Y coordinate of the shape's top left corner, relative to the origin of its parent object.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Y</em>' attribute.
     * @see #setY(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkShape_Y()
     * @model default="0.0" required="true"
     * @generated
     */
    double getY();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkShape#getY <em>Y</em>}' attribute.
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
     * Convenience method to set the shape's width and height simultaneously by calling their respective set methods.
     * <!-- end-model-doc -->
     * @model annotation="http://www.eclipse.org/emf/2002/GenModel body='setWidth(width);\nsetHeight(height);'"
     * @generated
     */
    void setDimensions(double width, double height);

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Convenience method to set the shape's x and y coordinates simultaneously by calling their respective set methods.
     * <!-- end-model-doc -->
     * @model annotation="http://www.eclipse.org/emf/2002/GenModel body='setX(x);\nsetY(y);'"
     * @generated
     */
    void setLocation(double x, double y);

} // ElkShape
