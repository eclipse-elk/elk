/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData;

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Md Group</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdGroup#getChildren <em>Children</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdGroup()
 * @model
 * @generated
 */
public interface MdGroup extends MdGroupOrOption
{
  /**
   * Returns the value of the '<em><b>Children</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdGroupOrOption}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Children</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Children</em>' containment reference list.
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdGroup_Children()
   * @model containment="true"
   * @generated
   */
  EList<MdGroupOrOption> getChildren();

} // MdGroup
