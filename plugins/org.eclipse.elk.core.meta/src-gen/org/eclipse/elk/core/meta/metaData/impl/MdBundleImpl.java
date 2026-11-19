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
package org.eclipse.elk.core.meta.metaData.impl;

import java.util.Collection;

import org.eclipse.elk.core.meta.metaData.MdBundle;
import org.eclipse.elk.core.meta.metaData.MdBundleMember;
import org.eclipse.elk.core.meta.metaData.MetaDataPackage;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

import org.eclipse.emf.ecore.util.EObjectContainmentEList;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Md Bundle</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl#getTargetClass <em>Target Class</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl#getDocumentationFolder <em>Documentation Folder</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl#getIdPrefix <em>Id Prefix</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl#getMembers <em>Members</em>}</li>
 * </ul>
 *
 * @generated
 */
public class MdBundleImpl extends MinimalEObjectImpl.Container implements MdBundle
{
  /**
   * The default value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected static final String LABEL_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected String label = LABEL_EDEFAULT;

  /**
   * The default value of the '{@link #getTargetClass() <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetClass()
   * @generated
   * @ordered
   */
  protected static final String TARGET_CLASS_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getTargetClass() <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetClass()
   * @generated
   * @ordered
   */
  protected String targetClass = TARGET_CLASS_EDEFAULT;

  /**
   * The default value of the '{@link #getDocumentationFolder() <em>Documentation Folder</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDocumentationFolder()
   * @generated
   * @ordered
   */
  protected static final String DOCUMENTATION_FOLDER_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getDocumentationFolder() <em>Documentation Folder</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDocumentationFolder()
   * @generated
   * @ordered
   */
  protected String documentationFolder = DOCUMENTATION_FOLDER_EDEFAULT;

  /**
   * The default value of the '{@link #getIdPrefix() <em>Id Prefix</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getIdPrefix()
   * @generated
   * @ordered
   */
  protected static final String ID_PREFIX_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getIdPrefix() <em>Id Prefix</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getIdPrefix()
   * @generated
   * @ordered
   */
  protected String idPrefix = ID_PREFIX_EDEFAULT;

  /**
   * The cached value of the '{@link #getMembers() <em>Members</em>}' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMembers()
   * @generated
   * @ordered
   */
  protected EList<MdBundleMember> members;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected MdBundleImpl()
  {
    super();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  protected EClass eStaticClass()
  {
    return MetaDataPackage.Literals.MD_BUNDLE;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getLabel()
  {
    return label;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLabel(String newLabel)
  {
    String oldLabel = label;
    label = newLabel;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_BUNDLE__LABEL, oldLabel, label));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getTargetClass()
  {
    return targetClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setTargetClass(String newTargetClass)
  {
    String oldTargetClass = targetClass;
    targetClass = newTargetClass;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_BUNDLE__TARGET_CLASS, oldTargetClass, targetClass));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getDocumentationFolder()
  {
    return documentationFolder;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDocumentationFolder(String newDocumentationFolder)
  {
    String oldDocumentationFolder = documentationFolder;
    documentationFolder = newDocumentationFolder;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_BUNDLE__DOCUMENTATION_FOLDER, oldDocumentationFolder, documentationFolder));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getIdPrefix()
  {
    return idPrefix;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setIdPrefix(String newIdPrefix)
  {
    String oldIdPrefix = idPrefix;
    idPrefix = newIdPrefix;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_BUNDLE__ID_PREFIX, oldIdPrefix, idPrefix));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<MdBundleMember> getMembers()
  {
    if (members == null)
    {
      members = new EObjectContainmentEList<MdBundleMember>(MdBundleMember.class, this, MetaDataPackage.MD_BUNDLE__MEMBERS);
    }
    return members;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public NotificationChain eInverseRemove(InternalEObject otherEnd, int featureID, NotificationChain msgs)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_BUNDLE__MEMBERS:
        return ((InternalEList<?>)getMembers()).basicRemove(otherEnd, msgs);
    }
    return super.eInverseRemove(otherEnd, featureID, msgs);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public Object eGet(int featureID, boolean resolve, boolean coreType)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_BUNDLE__LABEL:
        return getLabel();
      case MetaDataPackage.MD_BUNDLE__TARGET_CLASS:
        return getTargetClass();
      case MetaDataPackage.MD_BUNDLE__DOCUMENTATION_FOLDER:
        return getDocumentationFolder();
      case MetaDataPackage.MD_BUNDLE__ID_PREFIX:
        return getIdPrefix();
      case MetaDataPackage.MD_BUNDLE__MEMBERS:
        return getMembers();
    }
    return super.eGet(featureID, resolve, coreType);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @SuppressWarnings("unchecked")
  @Override
  public void eSet(int featureID, Object newValue)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_BUNDLE__LABEL:
        setLabel((String)newValue);
        return;
      case MetaDataPackage.MD_BUNDLE__TARGET_CLASS:
        setTargetClass((String)newValue);
        return;
      case MetaDataPackage.MD_BUNDLE__DOCUMENTATION_FOLDER:
        setDocumentationFolder((String)newValue);
        return;
      case MetaDataPackage.MD_BUNDLE__ID_PREFIX:
        setIdPrefix((String)newValue);
        return;
      case MetaDataPackage.MD_BUNDLE__MEMBERS:
        getMembers().clear();
        getMembers().addAll((Collection<? extends MdBundleMember>)newValue);
        return;
    }
    super.eSet(featureID, newValue);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public void eUnset(int featureID)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_BUNDLE__LABEL:
        setLabel(LABEL_EDEFAULT);
        return;
      case MetaDataPackage.MD_BUNDLE__TARGET_CLASS:
        setTargetClass(TARGET_CLASS_EDEFAULT);
        return;
      case MetaDataPackage.MD_BUNDLE__DOCUMENTATION_FOLDER:
        setDocumentationFolder(DOCUMENTATION_FOLDER_EDEFAULT);
        return;
      case MetaDataPackage.MD_BUNDLE__ID_PREFIX:
        setIdPrefix(ID_PREFIX_EDEFAULT);
        return;
      case MetaDataPackage.MD_BUNDLE__MEMBERS:
        getMembers().clear();
        return;
    }
    super.eUnset(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public boolean eIsSet(int featureID)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_BUNDLE__LABEL:
        return LABEL_EDEFAULT == null ? label != null : !LABEL_EDEFAULT.equals(label);
      case MetaDataPackage.MD_BUNDLE__TARGET_CLASS:
        return TARGET_CLASS_EDEFAULT == null ? targetClass != null : !TARGET_CLASS_EDEFAULT.equals(targetClass);
      case MetaDataPackage.MD_BUNDLE__DOCUMENTATION_FOLDER:
        return DOCUMENTATION_FOLDER_EDEFAULT == null ? documentationFolder != null : !DOCUMENTATION_FOLDER_EDEFAULT.equals(documentationFolder);
      case MetaDataPackage.MD_BUNDLE__ID_PREFIX:
        return ID_PREFIX_EDEFAULT == null ? idPrefix != null : !ID_PREFIX_EDEFAULT.equals(idPrefix);
      case MetaDataPackage.MD_BUNDLE__MEMBERS:
        return members != null && !members.isEmpty();
    }
    return super.eIsSet(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    if (eIsProxy()) return super.toString();

    StringBuilder result = new StringBuilder(super.toString());
    result.append(" (label: ");
    result.append(label);
    result.append(", targetClass: ");
    result.append(targetClass);
    result.append(", documentationFolder: ");
    result.append(documentationFolder);
    result.append(", idPrefix: ");
    result.append(idPrefix);
    result.append(')');
    return result.toString();
  }

} //MdBundleImpl
