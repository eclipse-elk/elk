/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Configuration</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSamples <em>Samples</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getForm <em>Form</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getNodes <em>Nodes</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getEdges <em>Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMW <em>MW</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxWidth <em>Max Width</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMD <em>MD</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxDegree <em>Max Degree</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isPF <em>PF</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFraction <em>Fraction</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getHierarchy <em>Hierarchy</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSeed <em>Seed</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFormat <em>Format</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFilename <em>Filename</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration()
 * @model
 * @generated
 */
public interface Configuration extends EObject
{
  /**
   * Returns the value of the '<em><b>Samples</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Samples</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Samples</em>' attribute.
   * @see #setSamples(int)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Samples()
   * @model
   * @generated
   */
  int getSamples();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSamples <em>Samples</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Samples</em>' attribute.
   * @see #getSamples()
   * @generated
   */
  void setSamples(int value);

  /**
   * Returns the value of the '<em><b>Form</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.core.debug.grandom.gRandom.Form}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Form</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Form</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Form
   * @see #setForm(Form)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Form()
   * @model
   * @generated
   */
  Form getForm();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getForm <em>Form</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Form</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Form
   * @see #getForm()
   * @generated
   */
  void setForm(Form value);

  /**
   * Returns the value of the '<em><b>Nodes</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Nodes</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Nodes</em>' containment reference.
   * @see #setNodes(Nodes)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Nodes()
   * @model containment="true"
   * @generated
   */
  Nodes getNodes();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getNodes <em>Nodes</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Nodes</em>' containment reference.
   * @see #getNodes()
   * @generated
   */
  void setNodes(Nodes value);

  /**
   * Returns the value of the '<em><b>Edges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Edges</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Edges</em>' containment reference.
   * @see #setEdges(Edges)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Edges()
   * @model containment="true"
   * @generated
   */
  Edges getEdges();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getEdges <em>Edges</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Edges</em>' containment reference.
   * @see #getEdges()
   * @generated
   */
  void setEdges(Edges value);

  /**
   * Returns the value of the '<em><b>MW</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>MW</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>MW</em>' attribute.
   * @see #setMW(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_MW()
   * @model
   * @generated
   */
  boolean isMW();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMW <em>MW</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>MW</em>' attribute.
   * @see #isMW()
   * @generated
   */
  void setMW(boolean value);

  /**
   * Returns the value of the '<em><b>Max Width</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Max Width</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Max Width</em>' attribute.
   * @see #setMaxWidth(Integer)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_MaxWidth()
   * @model
   * @generated
   */
  Integer getMaxWidth();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxWidth <em>Max Width</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Max Width</em>' attribute.
   * @see #getMaxWidth()
   * @generated
   */
  void setMaxWidth(Integer value);

  /**
   * Returns the value of the '<em><b>MD</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>MD</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>MD</em>' attribute.
   * @see #setMD(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_MD()
   * @model
   * @generated
   */
  boolean isMD();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isMD <em>MD</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>MD</em>' attribute.
   * @see #isMD()
   * @generated
   */
  void setMD(boolean value);

  /**
   * Returns the value of the '<em><b>Max Degree</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Max Degree</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Max Degree</em>' attribute.
   * @see #setMaxDegree(Integer)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_MaxDegree()
   * @model
   * @generated
   */
  Integer getMaxDegree();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getMaxDegree <em>Max Degree</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Max Degree</em>' attribute.
   * @see #getMaxDegree()
   * @generated
   */
  void setMaxDegree(Integer value);

  /**
   * Returns the value of the '<em><b>PF</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>PF</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>PF</em>' attribute.
   * @see #setPF(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_PF()
   * @model
   * @generated
   */
  boolean isPF();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#isPF <em>PF</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>PF</em>' attribute.
   * @see #isPF()
   * @generated
   */
  void setPF(boolean value);

  /**
   * Returns the value of the '<em><b>Fraction</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Fraction</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Fraction</em>' containment reference.
   * @see #setFraction(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Fraction()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getFraction();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFraction <em>Fraction</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Fraction</em>' containment reference.
   * @see #getFraction()
   * @generated
   */
  void setFraction(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Hierarchy</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Hierarchy</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Hierarchy</em>' containment reference.
   * @see #setHierarchy(Hierarchy)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Hierarchy()
   * @model containment="true"
   * @generated
   */
  Hierarchy getHierarchy();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getHierarchy <em>Hierarchy</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Hierarchy</em>' containment reference.
   * @see #getHierarchy()
   * @generated
   */
  void setHierarchy(Hierarchy value);

  /**
   * Returns the value of the '<em><b>Seed</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Seed</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Seed</em>' attribute.
   * @see #setSeed(Integer)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Seed()
   * @model
   * @generated
   */
  Integer getSeed();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getSeed <em>Seed</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Seed</em>' attribute.
   * @see #getSeed()
   * @generated
   */
  void setSeed(Integer value);

  /**
   * Returns the value of the '<em><b>Format</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.core.debug.grandom.gRandom.Formats}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Format</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Format</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Formats
   * @see #setFormat(Formats)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Format()
   * @model
   * @generated
   */
  Formats getFormat();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFormat <em>Format</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Format</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Formats
   * @see #getFormat()
   * @generated
   */
  void setFormat(Formats value);

  /**
   * Returns the value of the '<em><b>Filename</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Filename</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Filename</em>' attribute.
   * @see #setFilename(String)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getConfiguration_Filename()
   * @model
   * @generated
   */
  String getFilename();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration#getFilename <em>Filename</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Filename</em>' attribute.
   * @see #getFilename()
   * @generated
   */
  void setFilename(String value);

} // Configuration
