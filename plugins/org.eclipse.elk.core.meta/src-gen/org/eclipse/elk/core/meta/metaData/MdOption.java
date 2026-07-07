/**
 * ******************************************************************************
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
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData;

import org.eclipse.emf.common.util.EList;

import org.eclipse.xtext.common.types.JvmTypeReference;

import org.eclipse.xtext.xbase.XExpression;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Md Option</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#isDeprecated <em>Deprecated</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#isAdvanced <em>Advanced</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#isProgrammatic <em>Programmatic</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#isOutput <em>Output</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#isGlobal <em>Global</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getType <em>Type</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getDescription <em>Description</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getDefaultValue <em>Default Value</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getLowerBound <em>Lower Bound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getUpperBound <em>Upper Bound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getTargets <em>Targets</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getLegacyIds <em>Legacy Ids</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOption#getDependencies <em>Dependencies</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption()
 * @model
 * @generated
 */
public interface MdOption extends MdGroupOrOption
{
  /**
   * Returns the value of the '<em><b>Deprecated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Deprecated</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Deprecated</em>' attribute.
   * @see #setDeprecated(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Deprecated()
   * @model
   * @generated
   */
  boolean isDeprecated();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#isDeprecated <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Deprecated</em>' attribute.
   * @see #isDeprecated()
   * @generated
   */
  void setDeprecated(boolean value);

  /**
   * Returns the value of the '<em><b>Advanced</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Advanced</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Advanced</em>' attribute.
   * @see #setAdvanced(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Advanced()
   * @model
   * @generated
   */
  boolean isAdvanced();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#isAdvanced <em>Advanced</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Advanced</em>' attribute.
   * @see #isAdvanced()
   * @generated
   */
  void setAdvanced(boolean value);

  /**
   * Returns the value of the '<em><b>Programmatic</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Programmatic</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Programmatic</em>' attribute.
   * @see #setProgrammatic(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Programmatic()
   * @model
   * @generated
   */
  boolean isProgrammatic();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#isProgrammatic <em>Programmatic</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Programmatic</em>' attribute.
   * @see #isProgrammatic()
   * @generated
   */
  void setProgrammatic(boolean value);

  /**
   * Returns the value of the '<em><b>Output</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Output</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Output</em>' attribute.
   * @see #setOutput(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Output()
   * @model
   * @generated
   */
  boolean isOutput();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#isOutput <em>Output</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Output</em>' attribute.
   * @see #isOutput()
   * @generated
   */
  void setOutput(boolean value);

  /**
   * Returns the value of the '<em><b>Global</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Global</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Global</em>' attribute.
   * @see #setGlobal(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Global()
   * @model
   * @generated
   */
  boolean isGlobal();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#isGlobal <em>Global</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Global</em>' attribute.
   * @see #isGlobal()
   * @generated
   */
  void setGlobal(boolean value);

  /**
   * Returns the value of the '<em><b>Type</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Type</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Type</em>' containment reference.
   * @see #setType(JvmTypeReference)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Type()
   * @model containment="true"
   * @generated
   */
  JvmTypeReference getType();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getType <em>Type</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Type</em>' containment reference.
   * @see #getType()
   * @generated
   */
  void setType(JvmTypeReference value);

  /**
   * Returns the value of the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Label</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Label</em>' attribute.
   * @see #setLabel(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Label()
   * @model
   * @generated
   */
  String getLabel();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getLabel <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Label</em>' attribute.
   * @see #getLabel()
   * @generated
   */
  void setLabel(String value);

  /**
   * Returns the value of the '<em><b>Description</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Description</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Description</em>' attribute.
   * @see #setDescription(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Description()
   * @model
   * @generated
   */
  String getDescription();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getDescription <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Description</em>' attribute.
   * @see #getDescription()
   * @generated
   */
  void setDescription(String value);

  /**
   * Returns the value of the '<em><b>Default Value</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Default Value</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Default Value</em>' containment reference.
   * @see #setDefaultValue(XExpression)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_DefaultValue()
   * @model containment="true"
   * @generated
   */
  XExpression getDefaultValue();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getDefaultValue <em>Default Value</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Default Value</em>' containment reference.
   * @see #getDefaultValue()
   * @generated
   */
  void setDefaultValue(XExpression value);

  /**
   * Returns the value of the '<em><b>Lower Bound</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Lower Bound</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Lower Bound</em>' containment reference.
   * @see #setLowerBound(XExpression)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_LowerBound()
   * @model containment="true"
   * @generated
   */
  XExpression getLowerBound();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getLowerBound <em>Lower Bound</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Lower Bound</em>' containment reference.
   * @see #getLowerBound()
   * @generated
   */
  void setLowerBound(XExpression value);

  /**
   * Returns the value of the '<em><b>Upper Bound</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Upper Bound</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Upper Bound</em>' containment reference.
   * @see #setUpperBound(XExpression)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_UpperBound()
   * @model containment="true"
   * @generated
   */
  XExpression getUpperBound();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOption#getUpperBound <em>Upper Bound</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Upper Bound</em>' containment reference.
   * @see #getUpperBound()
   * @generated
   */
  void setUpperBound(XExpression value);

  /**
   * Returns the value of the '<em><b>Targets</b></em>' attribute list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdOptionTargetType}.
   * The literals are from the enumeration {@link org.eclipse.elk.core.meta.metaData.MdOptionTargetType}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Targets</em>' attribute list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Targets</em>' attribute list.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionTargetType
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Targets()
   * @model unique="false"
   * @generated
   */
  EList<MdOptionTargetType> getTargets();

  /**
   * Returns the value of the '<em><b>Legacy Ids</b></em>' attribute list.
   * The list contents are of type {@link java.lang.String}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Legacy Ids</em>' attribute list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Legacy Ids</em>' attribute list.
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_LegacyIds()
   * @model unique="false"
   * @generated
   */
  EList<String> getLegacyIds();

  /**
   * Returns the value of the '<em><b>Dependencies</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdOptionDependency}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Dependencies</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Dependencies</em>' containment reference list.
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOption_Dependencies()
   * @model containment="true"
   * @generated
   */
  EList<MdOptionDependency> getDependencies();

} // MdOption
